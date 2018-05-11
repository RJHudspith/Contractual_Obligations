/**
   @file nrqcd.c
   @brief start of our "on the fly" NRQCD code

   TODO :: Really important that we think about multiple 
   time sources as I think we might be in real trouble with this
   and probably will have to compute the whole propagator and 
   put into memory
 **/
#include "common.h"

#include "clover.h"
#include "GLU_timer.h"
#include "halfspinor_ops.h"
#include "matrix_ops.h"
#include "progress_bar.h"
#include "sources.h"
#include "spin_dependent.h"
#include "spin_independent.h"

// so in principle this could be checkerboarded and the
// whole LCU loop could be done in this level in parallel
static void
evolve_H( struct NRQCD_fields *F ,
	  const size_t t ,
	  const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C0 ) < 1E-12 ) return ;

  const double C0 = -NRQCD.C0 / ( 2. * NRQCD.M_0 ) ;

  // simplest hamiltonian term first \grad^2 / 2M_0
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {

    // temporary storage matrices
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t idx = i + t*LCU ;
    size_t mu , d ;
    
    // set hamiltonian storage to zero
    zero_colormatrix( F -> H[i].D[0] ) ;
    zero_colormatrix( F -> H[i].D[1] ) ;
    zero_colormatrix( F -> H[i].D[2] ) ;
    zero_colormatrix( F -> H[i].D[3] ) ;

    // inner mu sum much more cache friendly
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t fwd  = lat[ i ].neighbor[mu] ;
      const size_t bck  = lat[ i ].back[mu] ;
      const size_t bck2 = lat[idx].back[mu] ;
      
      for( d = 0 ; d < NS ; d++ ) {
	multab( (void*)A , (void*)lat[idx].O[mu] , (void*)F->S[fwd].D[d] ) ;
	multabdag( (void*)B , (void*)lat[bck2].O[mu] , (void*)F->S[bck].D[d] ) ;
	// A = A + B 
	add_mat( (void*)A , (void*)B ) ;
	// A = A + B - 2 F -> S[i].D[d]
	colormatrix_Saxpy( A , F->S[i].D[d], -2. ) ;
	// F = C0 * ( A + B - 2 F -> S[i].D[d] )
	colormatrix_Saxpy( F -> H[i].D[d] , A , C0 ) ;
      }
    }
  }

  halfspinor_Saxpy( F -> S , F -> H , -1./(2.*NRQCD.N) ) ;
    
  return ;
}

// applies the hamiltonian s.t. S = ( 1 - H/2 ) S
static void
evolve_dH( struct NRQCD_fields *F ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD )
{
  zero_halfspinor( F -> H ) ;
  
  // atomically accumulate result into F -> H
  term_C1_C6( F , t , NRQCD ) ;
  term_C2( F , t , NRQCD ) ;
  term_C3( F , t , NRQCD ) ;
  term_C4( F , t , NRQCD ) ;
  term_C5( F , t , NRQCD ) ;
  term_C7( F , t , NRQCD ) ;
  term_C8( F , t , NRQCD ) ;
  term_C9EB( F , t , NRQCD ) ;
  term_C10EB( F , t , NRQCD ) ;
  term_C11( F , t , NRQCD ) ;
  
  halfspinor_Saxpy( F -> S , F -> H , -1. ) ;
    
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

  evolve_dH( F , t , NRQCD ) ;

  // evolve just with C0 term
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD ) ;
  }
  // mutliply by temporal link
  for( i = 0 ; i < LCU ; i++ ) {    
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    const size_t idx = i + t*LCU ;
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      colormatrix_equiv( (void*)A , (void*)F -> S[i].D[d] ) ;
      multabdag( (void*)F -> S[i].D[d] , (void*)lat[idx].O[ ND-1 ] , (void*)A ) ; 
    }
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
  for( i = 0 ; i < LCU ; i++ ) {    
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    const size_t idx = lat[i+t*LCU].back[ ND-1 ] ;
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      colormatrix_equiv( (void*)A , (void*)F -> S[i].D[d] ) ;
      multab( (void*)F -> S[i].D[d] , (void*)lat[idx].O[ ND-1 ] , (void*)A ) ; 
    }
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
  size_t i , mu , j ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[i].O[mu][j] *= tadpole ;
      }
    }
  }
  return ;
}

// compute our various NRQCD propagators
int
compute_nrqcd_props( struct propagator *prop ,
		     const size_t nprops )
{
  // loop propagators checking to see if any are non-rel
  // check tadpole factors are the same
  double tadref = 0.0 ;
  GLU_bool FLY_NREL = GLU_FALSE ;
  size_t n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].basis == NREL_CORR ) {
      FLY_NREL = GLU_TRUE ;
      if( fabs( tadref ) < 1E-12 ) {
	tadref = prop[n].NRQCD.U0 ;
      } else {
	if( fabs( prop[n].NRQCD.U0 - tadref ) > 1E-12 ||
	    prop[n].NRQCD.U0 < 1E-12 ) {
	  fprintf( stderr , "[NRQCD] poor selection of tadpole factor "
		   "prop %zu || fac %f \n" , n , prop[n].NRQCD.U0 ) ;
	  return FAILURE ;
	}
      }
    }
  }
  // if we aren't doing any NRQCD props then we successfully do nothing
  if( FLY_NREL == GLU_FALSE ) {
    fprintf( stdout , "[NRQCD] Not computing NRQCD props on the fly\n" ) ;
    return SUCCESS ;
  }
  // check for an empty gauge field
  if( lat == NULL ) {
    fprintf( stderr , "[NRQCD] Empty gauge field, add with -c on command line\n" ) ;
    return FAILURE ;
  }  
  // otherwise we initialise all this gubbins
  struct NRQCD_fields F ;
  
  corr_malloc( (void**)&F.S  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.S1 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.S2 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.S3 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.S4 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.H  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) ;
  corr_malloc( (void**)&F.Fmunu , ALIGNMENT , LCU*sizeof( struct field ) ) ;
  
  size_t i , j ;
  for( i = 0 ; i < LCU ; i++ ) {
    corr_malloc( (void**)&F.Fmunu[i].O , ALIGNMENT , 6*sizeof( double complex* ) ) ;
    for( j = 0 ; j < 6 ; j++ ) {
      corr_malloc( (void**)&F.Fmunu[i].O[j] , ALIGNMENT , NCNC*sizeof( double complex ) ) ;
    }
  }

  // initialise the timer
  start_timer() ;

  // do the tadpole improvement on the gauge field
  tadpole_improve( 1./tadref ) ;

  // loop N props this far out as we might want to have different source
  // positions
  for( n = 0 ; n < nprops ; n++ ) {

    if( prop[n].basis != NREL_CORR ) continue ;
    
    // allocate the heavy propagator
    corr_malloc( (void**)&prop[n].H , ALIGNMENT , LVOLUME*sizeof( struct halfspinor ) ) ;
    
    // set up the source
    initialise_source( prop[n].H , prop[n].source , prop[n].origin ) ;
    
    // evolve for all timeslices can we parallelise this? - nope
    size_t t , tnew , tprev = ( prop[n].origin[ND-1] )%LT ;
    for( t = 0 ; t < LT-1 ; t++ ) {
      
      // are we going backward or forward?
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	tnew = ( tprev + LT - 1 )%LT ;
      } else {
	tnew = ( tprev + LT + 1 )%LT ;
      }

      // compute the lattice clovers
      compute_clovers( &F , lat , tprev ) ;

      // do a copy in here
      for( i = 0 ; i < LCU ; i++ ) {
	const size_t idx = LCU*tprev ;
	size_t k ;
	for( k = 0 ; k < NCNC ; k++ ) {
	  F.S[i].D[0][k] = prop[n].H[i+idx].D[0][k] ;
	  F.S[i].D[1][k] = prop[n].H[i+idx].D[1][k] ;
	  F.S[i].D[2][k] = prop[n].H[i+idx].D[2][k] ;
	  F.S[i].D[3][k] = prop[n].H[i+idx].D[3][k] ;
	}
      }
      
      // compute the propagator
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	nrqcd_prop_bwd( &F , tprev , prop[n].NRQCD ) ;
      } else {
	nrqcd_prop_fwd( &F , tprev , prop[n].NRQCD ) ;
      }
      
      // set G(t+1) from temp1
      for( i = 0 ; i < LCU ; i++ ) {
	const size_t shift = LCU*tnew ;
	size_t k ;
	for( k = 0 ; k < NCNC ; k++ ) {
	  prop[n].H[i+shift].D[0][k] = F.S[i].D[0][k] ;
	  prop[n].H[i+shift].D[1][k] = F.S[i].D[1][k] ;
	  prop[n].H[i+shift].D[2][k] = F.S[i].D[2][k] ;
	  prop[n].H[i+shift].D[3][k] = F.S[i].D[3][k] ;
	}
      }

      // set the time index
      tprev = tnew ;

      progress_bar( t , LT-1 ) ;
    }
    // end the loop on propagators
  }

  // tadpole unimprove the gauge field
  tadpole_improve( tadref ) ;
  
  // tell us the time
  print_time() ;

  // memory frees
  free( F.S ) ; free( F.H ) ; free( F.S1 ) ;
  free( F.S2 ) ; free( F.S3 ) ; free( F.S4 ) ;

  // free the field strength tensor
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < 6 ; j++ ) {
      free( F.Fmunu[i].O[j] ) ;
    }
    free( F.Fmunu[i].O ) ;
  }
  free( F.Fmunu ) ;
  
  return SUCCESS ;
}

int
free_nrqcd_props( struct propagator *prop ,
		  const size_t nprops )
{
  size_t n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].H != NULL ) {
      free( prop[n].H ) ;
    }
  }
  return SUCCESS ;
}
