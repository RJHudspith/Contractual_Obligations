/**
   @file evolve.c
   @brief evolve the NRQCD action
 */
#include "common.h"

#include "clover.h"           // O(a^2) improved clover
#include "GLU_timer.h"        // tells us how long it takes
#include "halfspinor_ops.h"   // halfspinor_Saxpy
#include "matrix_ops.h"       // colormatrix_*
#include "mmul.h"             // multabs
#include "progress_bar.h"     // show the progress bar
#include "sources.h"          // initialise the source
#include "spin_dependent.h"   // spin-dependent NRQCD terms
#include "spin_independent.h" // spin-independent NRQCD terms

// so in principle this could be checkerboarded and the
// whole LCU loop could be done in this level in parallel
static void
evolve_H( struct NRQCD_fields *F ,
	  const size_t t ,
	  const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C0 ) < NRQCD_TOL ) return ;

  const double C0 = -NRQCD.C0 / ( 2. * NRQCD.M_0 ) ;

  // simplest hamiltonian term first \grad^2 / 2M_0
  // grad^2 is inlined from derivs.c
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // temporary storage matrices
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t Uidx = i + t*LCU ;
    size_t mu , d ;
    
    // set hamiltonian storage to zero
    zero_halfspinor( &F -> H[i] ) ;

    // inner mu sum much more cache friendly
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      
      const size_t Sfwd = lat[ i ].neighbor[mu] ;
      const size_t Sbck = lat[ i ].back[mu] ;
      const size_t Ubck = lat[ Uidx ].back[mu] ;

      for( d = 0 ; d < NS ; d++ ) {
	// computes A = U(x) S(x+\mu)
	multab( (void*)A ,
		(void*)lat[ Uidx ].O[mu] ,
		(void*)F -> S[ Sfwd ].D[d] ) ;
	// computes B = U^\dag(x-\mu) S(x-\mu)
	multabdag( (void*)B ,
		   (void*)lat[ Ubck ].O[mu] ,
		   (void*)F -> S[ Sbck ].D[d] ) ;
	// A = A + B 
	add_mat( (void*)A , (void*)B ) ;
	// A = A - 2 F -> S[i].D[d]
	colormatrix_Saxpy( A , F -> S[i].D[d] , -2. ) ;
	// der2 = A + der2
	colormatrix_Saxpy( F -> H[i].D[d] , A , C0 ) ;
      }
    }
  }

#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    halfspinor_Saxpy( &F -> S[i] , F -> H[i] , -1./(2.*NRQCD.N) ) ;
  }
  return ;
}

// applies the hamiltonian s.t. S = ( 1 - dH ) S
static void
evolve_dH( struct NRQCD_fields *F ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    zero_halfspinor( &F -> H[i] ) ;
  
    // atomically accumulate result into F -> H
    term_C1_C6( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C2( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C3( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C4( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C5( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C9EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C10EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;
  }

  // these three are really tricky to loop-fuse into the above
  term_C7( F , t , NRQCD ) ;
  term_C8( F , t , NRQCD ) ;
  term_C11( F , t , NRQCD ) ;

#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    halfspinor_Saxpy( &F -> S[i] , F -> H[i] , -1. ) ;
  }

  return ;
}

// writes out the result to S which is a LCU halfspinor
static int
nrqcd_prop_fwd( struct NRQCD_fields *F ,
		const size_t t ,
		const struct NRQCD_params NRQCD )
{  
  // loop power of NRQCD.N we apply the hamiltonian
  size_t n , i ;

  // only evolve the spin-dependent terms a little bit
  evolve_dH( F , t , NRQCD ) ;
  
  // evolve just with C0 term
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD ) ;
  }

  // mutliply by temporal link and put in temporary S1
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {    
    const size_t idx = i + t*LCU ;
    struct halfspinor res ;
    colormatrixdag_halfspinor( &res , lat[idx].O[ ND-1 ] , F -> S[i] ) ;
    F -> S[i] = res ;
  }

  // evolve again
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD ) ;
  }
  
  return SUCCESS ;
}

// writes out the result to S which is a LCU halfspinor
static int
nrqcd_prop_bwd( struct NRQCD_fields *F ,
		const size_t t ,
		const struct NRQCD_params NRQCD )
{  
  // loop power of NRQCD.N we apply the hamiltonian
  size_t n , i ;

  evolve_dH( F , t , NRQCD ) ;
  
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD ) ;
  }
    
  // mutliply by temporal link
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {    
    const size_t idx = lat[i+t*LCU].back[ ND-1 ] ;
    struct halfspinor res ;
    colormatrix_halfspinor( &res , lat[idx].O[ ND-1 ] , F -> S[i] ) ;
    F -> S[i] = res ;
  }

  // evolve again
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD ) ;
  }
  
  return SUCCESS ;
}

// tadpole improvement of the gauge field
static void
tadpole_improve( const double tadpole )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu , j ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[i].O[mu][j] *= tadpole ;
      }
    }
  }
  return ;
}

// this is the brains of the operation, evolves Hamiltonian from the
// source position
void
compute_props( struct propagator *prop ,
	       struct NRQCD_fields *F ,
	       const struct site *lat ,
	       const size_t nprops ,
	       const double tadref )
{
  // do the tadpole improvement on the gauge field
  tadpole_improve( 1./tadref ) ;

  size_t i , n ;
  // loop N props this far out as we might want to have different source
  // positions
  for( n = 0 ; n < nprops ; n++ ) {

    if( prop[n].basis != NREL_CORR ) continue ;

    size_t t , tnew , tprev = ( prop[n].origin[ND-1] )%LT ;

    // set up the source into F -> S
    initialise_source( F -> S ,
		       prop[n].source ,
		       prop[n].twists ,
		       prop[n].origin ) ;

    // do a copy in here
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      const size_t idx = LCU*tprev ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[0] , F -> S[i].D[0] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[1] , F -> S[i].D[1] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[2] , F -> S[i].D[2] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[3] , F -> S[i].D[3] ) ;
    }
    
    // evolve for all timeslices can we parallelise this? - nope
    for( t = 0 ; t < LT-1 ; t++ ) {
	
      // are we going backward or forward?
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	tnew = ( tprev + LT - 1 )%LT ;
      } else {
	tnew = ( tprev + LT + 1 )%LT ;
      }

      // compute the lattice clovers
      compute_clovers( F , lat , tprev , tadref ) ;
      	
      // compute the propagator
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	nrqcd_prop_bwd( F , tprev , prop[n].NRQCD ) ;
      } else {
	nrqcd_prop_fwd( F , tprev , prop[n].NRQCD ) ;
      }
    	
      // set G(t+1) from temp1
      #pragma omp for private(i)
      for( i = 0 ; i < LCU ; i++ ) {
	const size_t idx = LCU*tnew ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[0] , F -> S[i].D[0] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[1] , F -> S[i].D[1] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[2] , F -> S[i].D[2] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[3] , F -> S[i].D[3] ) ;
      }

      // set the time index
      tprev = tnew ;

      #pragma omp single nowait
      {
	progress_bar( t , LT-1 ) ;
      }
    }
  }
    
  // tadpole unimprove the gauge field
  tadpole_improve( tadref ) ;
 
  return ;
}
