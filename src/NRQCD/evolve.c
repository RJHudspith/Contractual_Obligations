/**
   @file evolve.c
   @brief evolve the NRQCD action
 */
#define _GNU_SOURCE // sincos

#include "common.h"

#include "clover.h"           // O(a^2) improved clover
#include "geometry.h"
#include "GLU_timer.h"        // tells us how long it takes
#include "halfspinor_ops.h"   // halfspinor_Saxpy
#include "matrix_ops.h"       // colormatrix_*
#include "mmul.h"             // multabs
#include "progress_bar.h"     // show the progress bar
#include "sources.h"          // initialise the source
#include "spin_dependent.h"   // spin-dependent NRQCD terms
#include "spin_independent.h" // spin-independent NRQCD terms

// for the tadpole improvement
#ifdef HAVE_IMMINTRIN_H
  #include "SSE2_OPS.h"

#if NC==3

// inline C0 function for SU3
#define inline_su3_C0()							\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ;\
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ;\
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;						\
  *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
  A++ ; B++ ; S++ ; H++ ;

#endif

// little inline computation of H = H + C0*( A + B - 2*S )
static inline void
C0_term( __m128d *H ,
	 const __m128d *A ,
	 const __m128d *B ,
	 const __m128d *S ,
	 const double C0 )
{
  register const __m128d m2 = _mm_setr_pd( -2 , -2 ) ;
  register const __m128d c0 = _mm_setr_pd( C0 , C0 ) ;
#if NC==3
  inline_su3_C0() ;
  inline_su3_C0() ;
  inline_su3_C0() ;
  inline_su3_C0() ;
#else
  size_t i ;
  for( i = 0 ; i < NS*NCNC ; i++ ) {
    *H = SSE2_FMA( c0 , SSE2_FMA( m2 , *S , _mm_add_pd( *A , *B ) ) , *H ) ; \
    A++ ; B++ ; S++ ; H++ ;
  }
#endif
  return ;
}

#else

// little inline computation of H = H + C0*( A + B - 2*S )
static inline void
C0_term( double complex *H ,
	 const double complex *A ,
	 const double complex *B ,
	 const double complex *S ,
	 const double C0 )
{
  size_t i ;
  for( i = 0 ; i < NS*NCNC ; i++ ) {
    *H += C0*( -2*(*S) + ( *A + *B ) ) ;
    A++ ; B++ ; S++ ; H++ ;
  }
  return ;
}

#endif

// inline mu loop macro for unrolling
#define inline__mu_C0(mu)					\
  Sfwd = lat[ i ].neighbor[mu] ;				\
  Sbck = lat[ i ].back[mu] ;					\
  Ubck = lat[ Uidx ].back[mu] ;					\
  dagger_gauge( C , (void*)lat[ Ubck ].O[mu] ) ;		\
  colormatrix_halfspinor( (void*)A.D ,				\
			  (const void*)lat[ Uidx ].O[mu] ,	\
			  (const void*)F -> S[ Sfwd ].D ) ;	\
  colormatrix_halfspinor( (void*)B.D , C ,			\
			  (const void*)F -> S[ Sbck ].D ) ;	\
  C0_term( (void*)F ->H[i].D , (void*)A.D , (void*)B.D ,	\
	   (void*)F ->S[i].D , C0 ) ;				\

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
    struct halfspinor A , B ;
    #ifdef HAVE_IMMINTRIN_H
    __m128d C[ NCNC ] ;
    #else
    double complex C[ NCNC ] ;
    #endif
    
    const size_t Uidx = i + t*LCU ;
    size_t Sfwd , Sbck , Ubck ;
    
    // set hamiltonian storage to zero
    zero_halfspinor( &F -> H[i] ) ;

    // inner mu sum much more cache friendly
    #if ND==4
    inline__mu_C0(0) ;
    inline__mu_C0(1) ;
    inline__mu_C0(2) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      inline__mu_C0(mu) ;
    }
    #endif

    // set S1 to S
    F -> S1[i] = F -> S[i] ;
    // S1 = S1 - H/(2n)
    halfspinor_Saxpy( &F -> S1[i] , F -> H[i] , -1./(2*NRQCD.N) ) ;
  }

  // shallow pointer swap
#pragma omp single
  {
    struct halfspinor *P = F -> S ;
    F -> S = F -> S1 ;
    F -> S1 = P ;
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

    term_C7( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C9EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C10EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;
  }

  // these three are really tricky to loop-fuse into the above
  term_C8( F , t , NRQCD ) ;
  term_C11( F , t , NRQCD ) ;

#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    #ifdef NRQCD_NONSYM
    halfspinor_Saxpy( &F -> S[i] , F -> H[i] , -1. ) ;
    #else
    halfspinor_Saxpy( &F -> S[i] , F -> H[i] , -0.5 ) ;
    #endif
  }

  return ;
}

// applies the hamiltonian s.t. S = ( 1 - dH ) S
static void
evolve_dH_fused( struct NRQCD_fields *F ,
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

    term_C7( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C9EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    term_C10EB( &F -> H[i] , F -> S , F -> Fmunu , i , t , NRQCD ) ;

    // set S1 to S
    F -> S1[i] = F -> S[i] ;
    #ifdef NRQCD_NONSYM
    halfspinor_Saxpy( &F -> S1[i] , F -> H[i] , -1. ) ;
    #else
    halfspinor_Saxpy( &F -> S1[i] , F -> H[i] , -0.5 ) ;
    #endif
  }
  
  // shallow pointer swap
#pragma omp single
  {
    struct halfspinor *P = F -> S ;
    F -> S = F -> S1 ;
    F -> S1 = P ;
  }
  
  return ;
}

// writes out the result to S which is a LCU halfspinor
static int
nrqcd_prop_fwd( struct NRQCD_fields *F ,
		const size_t t ,
		const struct NRQCD_params NRQCD ,
		const GLU_bool fuse_dH )
{
  // forward time
  const size_t tfwd = (t+LT+1)%Latt.dims[ND-1] ;
  
  // loop power of NRQCD.N we apply the hamiltonian
  size_t n , i ;

  // only evolve the spin-dependent terms a little bit
  if( fuse_dH == GLU_TRUE ) {
    evolve_dH_fused( F , t , NRQCD ) ;
  } else {
    evolve_dH( F , t , NRQCD ) ;
  }
  
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

  // update the clovers
  compute_clovers( F , lat , tfwd , NRQCD.U0 ) ;
  
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , tfwd , NRQCD ) ;
  }

#ifndef NRQCD_NONSYM
  // finally evolve the spin-dependent terms a little bit
  if( fuse_dH == GLU_TRUE ) {
    evolve_dH_fused( F , tfwd , NRQCD ) ;
  } else {
    evolve_dH( F , tfwd , NRQCD ) ;
  }
#endif
  
  return SUCCESS ;
}

// writes out the result to S which is a LCU halfspinor
static int
nrqcd_prop_bwd( struct NRQCD_fields *F ,
		const size_t t ,
		const struct NRQCD_params NRQCD ,
		const GLU_bool fuse_dH )
{  
  // loop power of NRQCD.N we apply the hamiltonian
  size_t n , i ;

  // when we go backwards the E-field gets flipped.
  // can change this just by the parameters C_2
  // and C_3 as they are the only ones that depend on E
  struct NRQCD_params NRQCD_flipped = NRQCD ;
  NRQCD_flipped.C2 *= -1 ;
  NRQCD_flipped.C3 *= -1 ;

  // only evolve the spin-dependent terms a little bit
  if( fuse_dH == GLU_TRUE ) {
    evolve_dH_fused( F , t , NRQCD_flipped ) ;
  } else {
    evolve_dH( F , t , NRQCD_flipped ) ;
  }
  
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , t , NRQCD_flipped ) ;
  }
    
  // mutliply by temporal link
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {    
    const size_t idx = lat[i+t*LCU].back[ ND-1 ] ;
    struct halfspinor res ;
    colormatrix_halfspinor( (void*)res.D ,
			    (const void*)lat[idx].O[ ND-1 ] ,
			    (const void*)F -> S[i].D ) ;
    F -> S[i] = res ;
  }

  const size_t tbck = (t+LT-1)%LT ;
  
  // update the clovers
  compute_clovers( F , lat , tbck , NRQCD.U0 ) ;

  // evolve again
  for( n = 0 ; n < NRQCD.N ; n++ ) {
    evolve_H( F , tbck , NRQCD_flipped ) ;
  }

#ifndef NRQCD_NONSYM
  // finally evolve the spin-dependent terms a little bit
  if( fuse_dH == GLU_TRUE ) {
    evolve_dH_fused( F , tbck , NRQCD_flipped ) ;
  } else {
    evolve_dH( F , tbck , NRQCD_flipped ) ;
  }
#endif
  
  return SUCCESS ;
}

// tadpole improvement of the gauge field
static void
tadpole_improve( struct site *lat ,
		 const double tadpole )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    #ifdef HAVE_IMMINTRIN_H
    __m128d *pU = (__m128d*)lat[i].O ;
    const register __m128d tad = _mm_setr_pd( tadpole , tadpole ) ;
    size_t j ;
    for( j = 0 ; j < ND*NCNC ; j++ ) {
      *pU = _mm_mul_pd( *pU , tad ) ;
      pU++ ;
    }
    #else
    size_t j , mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
        lat[i].O[mu][j] *= tadpole ;
      }
    }
    #endif
  }
  return ;
}

// accumulate twists
static void
apply_twist( struct site *lat ,
	     double sum_twist[ ND ] ,
	     const double this_twist[ ND ] ,
	     const double mom_source[ ND ] )
{
  double delta[ ND ] ;

  // compute the additional twist needed to apply to gauge field 
  size_t mu , i ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    const double new_mom = this_twist[mu] ;
    delta[ mu ] = ( new_mom ) - sum_twist[ mu ] ;
    // sum = sum + ( this - sum ) == this
    sum_twist[ mu ] = new_mom ; 
  }

  // twist only in the spatial directions
  for( mu = 0 ; mu < ND-1 ; mu++ ) {

    // skip if it is zero as it just multiplies by 1
    if( fabs( delta[ mu ] ) < NRQCD_TOL ) continue ;

    // phase is exp( -I * theta_\mu * 2 * M_PI / Latt.dims[mu] )
    // the factor of 2 is to keep it in line with momentum def
    // this is kosher as the twist is arbitrary so it can just
    // reabsorb this factor
    const double dlat = TWOPI*delta[mu]/Latt.dims[mu] ;
    double s , c ;
    sincos( dlat , &s , &c ) ;
    
    // this might be a rare place where we can nowait the loop
    #pragma omp for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {

      #ifdef HAVE_IMMINTRIN_H
      __m128d *pU = (__m128d*)lat[i].O[mu] ;
      const register __m128d ph = _mm_setr_pd( c , s ) ;
        #if NC == 3
        *pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; pU++ ;
        *pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; pU++ ;
        *pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; pU++ ;
	*pU = SSE2_MUL( *pU , ph ) ; 
        #else
        size_t j ;
        for( j = 0 ; j < NCNC ; j++ ) {
          *pU = SSE2_MUL( *pU , ph ) ; pU++ ;
        }
        #endif
      #else
      const double complex phase = c + I*s ;
      size_t j ;
      for( j = 0 ; j < NCNC ; j++ ) {
        lat[i].O[mu][j] *= phase ;
      }
      #endif
    }
  }
  
  return ;
}

// this is the brains of the operation, evolves Hamiltonian from the
// source position
void
compute_props( struct propagator *prop ,
	       struct NRQCD_fields *F ,
	       struct site *lat ,
	       const size_t nprops ,
	       const double tadref )
{
  size_t i , n ;

  // do the tadpole improvement on the gauge field
  tadpole_improve( lat , 1./tadref ) ;

  // initialise the accumulated sum of twists we have applied
  double sum_twist[ ND ] ;
  for( i = 0 ; i < ND ; i++ ) {
    sum_twist[ i ] = 0.0 ;
  }
  
  // loop N props this far out as we might want to have different source
  // positions
  for( n = 0 ; n < nprops ; n++ ) {

    if( prop[n].basis != NREL_CORR ) continue ;

    size_t t , tnew , tprev = ( prop[n].origin[ND-1] )%LT ;

    // if we can we fuse the loops in dH evolution to avoid a
    // barrier
    GLU_bool fuse_dH = GLU_FALSE ;
    if( fabs( prop[n].NRQCD.C8 )  < NRQCD_TOL &&
	fabs( prop[n].NRQCD.C11 ) < NRQCD_TOL ) {
      fuse_dH = GLU_TRUE ;
    }

    // apply a twist to the gauge field
    apply_twist( lat , sum_twist , prop[n].twist , prop[n].mom_source ) ;

    // set up the source into F -> S
    initialise_source( F -> S , F -> S1 , prop[n] ) ;

    // do a copy in here
    #pragma omp for nowait private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      const size_t idx = LCU*tprev ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[0] , F -> S[i].D[0] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[1] , F -> S[i].D[1] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[2] , F -> S[i].D[2] ) ;
      colormatrix_equiv_d2f( prop[n].H[i+idx].D[3] , F -> S[i].D[3] ) ;
    }

    // compute the initial lattice clovers
    compute_clovers( F , lat , tprev , tadref ) ;
    
    // evolve for all timeslices can we parallelise this? - nope
    for( t = 0 ; t < T_NRQCD-1 ; t++ ) {
	
      // are we going backward or forward?
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	tnew = ( tprev + LT - 1 )%LT ;
      } else {
	tnew = ( tprev + LT + 1 )%LT ;
      }
      	
      // compute the propagator
      if( prop[n].NRQCD.backward == GLU_TRUE ) {
	nrqcd_prop_bwd( F , tprev%LT , prop[n].NRQCD , fuse_dH ) ;
      } else {
	nrqcd_prop_fwd( F , tprev%LT , prop[n].NRQCD , fuse_dH ) ;
      }
    	
      // set G(t+1) from temp1
      const size_t idx = LCU*((t + 1 + ( prop[n].origin[ND-1] )%LT )%(T_NRQCD) ) ;
      
      #pragma omp for nowait private(i) 
      for( i = 0 ; i < LCU ; i++ ) {
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[0] , F -> S[i].D[0] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[1] , F -> S[i].D[1] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[2] , F -> S[i].D[2] ) ;
	colormatrix_equiv_d2f( prop[n].H[i+idx].D[3] , F -> S[i].D[3] ) ;
      }

      // set the time index
      tprev = tnew ;

      #pragma omp single nowait
      {
	progress_bar( t , T_NRQCD-1 ) ;
      }
    }
  }

  // untwist the gauge field
  double zero_twist[ ND ] ;
  for( i = 0 ; i < ND ; i++ ) {
    zero_twist[ i ] = 0.0 ;
  }
  apply_twist( lat , sum_twist , zero_twist , zero_twist ) ;

  // un-tadpole improve
  tadpole_improve( lat , tadref ) ;
  
  return ;
}
