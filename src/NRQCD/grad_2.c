/**
   @file grad_sq.c
   @brief routines for the second order derivatives
 */
#include "common.h"

#include "matrix_ops.h"      // add_mat() etc
#include "mmul.h"            // multab() etc
#include "halfspinor_ops.h"  // zero_halfspinor() etc

#ifdef HAVE_IMMINTRIN_H

#define inner_unroll_su3()\
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
    *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
    *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \
  *H = _mm_add_pd( SSE2_FMA( *S , f1 , _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) , _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) ) , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ; \

static inline void
sq_inner( __m128d *H ,
	  const __m128d *A ,
	  const __m128d *B ,
	  const __m128d *C ,
	  const __m128d *D ,
	  const __m128d *S )
{
  const __m128d f1 = _mm_set_pd( -2.5 , -2.5 ) ;
  const __m128d f2 = _mm_set_pd( 4/3. , 4/3. ) ;
  const __m128d f3 = _mm_set_pd( -1/12. , -1/12. ) ;
#if NC==3
  inner_unroll_su3() ;
  inner_unroll_su3() ;
  inner_unroll_su3() ;
  inner_unroll_su3() ;
#else
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H = _mm_add_pd( SSE2_FMA( *S , f1 ,
			       _mm_add_pd( _mm_mul_pd( f2 , _mm_add_pd( *A , *B ) ) ,
					   _mm_mul_pd( f3 , _mm_add_pd( *C , *D ) ) ) )      
		    , *H ) ; H++ ; S++ ; A++ ; B++ ; C++ ; D++ ;
  }
#endif
  return ;
}

#else

static inline void
sq_inner( double complex *H ,
	  const double complex *A ,
	  const double complex *B ,
	  const double complex *C ,
	  const double complex *D ,
	  const double complex *S )
{
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H += -5/2.*(*S) + 4/3.*(*A + *B) - 1/12.*(*C + *D) ;
    H++ ; S++ ; A++ ; B++ ; C++ ; D++ ;
  }
  return ;
}

#endif


// similar to grad, sum over all mu
void
gradsq( struct halfspinor *der2 ,
	const struct halfspinor *S ,
	const size_t i ,
	const size_t t )
{      
  zero_halfspinor( der2 ) ;

  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

  const size_t Uidx = i + t*LCU ;

  size_t mu , d ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Ubck = lat[ Uidx ].back[mu] ;
      
    for( d = 0 ; d < NS ; d++ ) {
      // computes A = U(x) S(x+\mu)
      multab( (void*)A ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      // computes B = U^\dag(x-\mu) S(x-\mu)
      multabdag( (void*)B ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      // A = A + B 
      add_mat( (void*)A , (void*)B ) ;
      // A = A - 2 F -> S[i].D[d]
      colormatrix_Saxpy( A , S[i].D[d], -2. ) ;
      // der2 = A + der2
      add_mat( (void*)der2 -> D[d] , (void*)A ) ;
    }
  }
  return ;
}

// same as grad, sum over all mu
void
grad_sq_LCU( struct halfspinor *der2 ,
	     const struct halfspinor *S ,
	     const size_t t )
{
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    zero_halfspinor( &der2[i] ) ;

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t Uidx = i + t*LCU ;

    size_t mu , d ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t Sfwd = lat[ i ].neighbor[mu] ;
      const size_t Sbck = lat[ i ].back[mu] ;
      const size_t Ubck = lat[ Uidx ].back[mu] ;
      
      for( d = 0 ; d < NS ; d++ ) {
	// computes A = U(x) S(x+\mu)
	multab( (void*)A ,
		(void*)lat[ Uidx ].O[mu] ,
		(void*)S[ Sfwd ].D[d] ) ;
	// computes B = U^\dag(x-\mu) S(x-\mu)
	multabdag( (void*)B ,
		   (void*)lat[ Ubck ].O[mu] ,
		   (void*)S[ Sbck ].D[d] ) ;
	// A = A + B 
	add_mat( (void*)A , (void*)B ) ;
	// A = A - 2 F -> S[i].D[d]
	colormatrix_Saxpy( A , S[i].D[d], -2. ) ;
	// der2 = A + der2
	add_mat( (void*)der2[i].D[d] , (void*)A ) ;
      }
    }
  }
  return ;
}

// improved \grad^2
void
gradsq_imp( struct halfspinor *der ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t )
{
  struct halfspinor A , B , C , D ;
#ifdef HAVE_IMMINTRIN_H
  __m128d b[ NCNC ] ;
#else
  double complex b[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
#endif
  
  const size_t Uidx = i + t*LCU ;

  zero_halfspinor( der ) ;

  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
  
    const size_t Sfwd  = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck  = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    colormatrix_halfspinor( (void*)A.D ,(void*)lat[ Uidx ].O[mu] ,
			    (void*)S[ Sfwd ].D ) ;
    dagger_gauge( (void*)b , (void*)lat[ Ubck ].O[mu] ) ;
    colormatrix_halfspinor( (void*)B.D , (void*)b ,(void*)S[ Sbck ].D ) ;

    colormatrix_halfspinor( (void*)C.D , (void*)Fmunu[i].O[6+2*mu] ,
			    (void*)S[ Sfwd2 ].D ) ;
    dagger_gauge( (void*)b , (void*)Fmunu[i].O[7+2*mu] ) ;
    colormatrix_halfspinor( (void*)D.D , (void*)b , (void*)S[ Sbck2 ].D ) ;

    // do the sum der = -5/2 S[i] + 4/3(A+B) - 1/12(C+D) 
    sq_inner( (void*)der -> D , (void*)A.D , (void*)B.D ,
	      (void*)C.D , (void*)D.D , (void*)S[i].D ) ;
  }
  return ;
}

// improved \grad^2
void
gradsq_imp_sigmaB( struct halfspinor *der ,
		   const struct halfspinor *S ,
		   const struct field *Fmunu ,
		   const size_t i ,
		   const size_t t )
{
#ifdef HAVE_IMMMINTRIN_H
  __mm128d b[ NCNC ] ;
#else
  double complex b[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
#endif
  
  const size_t Uidx = i + t*LCU ;

  zero_halfspinor( der ) ;

  // eww, lots of stack allocations here
  struct halfspinor sigmaB_S , Stmp , A , B , C , D ;

  sigmaB_halfspinor( &sigmaB_S , Fmunu[i] , S[i] );
  
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
  
    const size_t Sfwd  = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck  = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    sigmaB_halfspinor( &Stmp , Fmunu[Sfwd]  , S[Sfwd] );
    colormatrix_halfspinor( (void*)A.D , (void*)lat[ Uidx ].O[mu] ,
			    (void*)Stmp.D ) ;
    
    sigmaB_halfspinor( &Stmp  , Fmunu[Sbck]  , S[Sbck] );
    dagger_gauge( (void*)b , (void*)lat[ Ubck ].O[mu] ) ;
    colormatrix_halfspinor( (void*)B.D , (void*)b , (void*)Stmp.D ) ;

    sigmaB_halfspinor( &Stmp , Fmunu[Sfwd2] , S[Sfwd2] );
    colormatrix_halfspinor( (void*)C.D , (void*)Fmunu[i].O[6+2*mu] ,
			    (void*)Stmp.D ) ;
    
    sigmaB_halfspinor( &Stmp , Fmunu[Sbck2] , S[Sbck2] );
    dagger_gauge( (void*)b , (void*)Fmunu[i].O[7+2*mu] ) ;
    colormatrix_halfspinor( (void*)D.D , (void*)b , (void*)Stmp.D ) ;

    // do the sum der = -5/2 S[i] + 4/3(A+B) - 1/12(C+D) 
    sq_inner( (void*)der -> D , (void*)A.D , (void*)B.D ,
	      (void*)C.D , (void*)D.D , (void*)sigmaB_S.D ) ;
  }
  return ;
}

// improved \grad^2
void
sigmaB_gradsq_imp( struct halfspinor *der ,
		   const struct halfspinor *S ,
		   const struct field *Fmunu ,
		   const size_t i ,
		   const size_t t )
{
#ifdef HAVE_IMMINTRIN_H
  __m128d b[ NCNC ] ;
#else
  double complex b[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
#endif
  const size_t Uidx = i + t*LCU ;

  struct halfspinor res , A , B , C , D ; 
  zero_halfspinor( &res ) ;

  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
  
    const size_t Sfwd  = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck  = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    colormatrix_halfspinor( (void*)A.D ,(void*)lat[ Uidx ].O[mu] ,
			    (void*)S[ Sfwd ].D ) ;
    dagger_gauge( (void*)b , (void*)lat[ Ubck ].O[mu] ) ;
    colormatrix_halfspinor( (void*)B.D , (void*)b ,(void*)S[ Sbck ].D ) ;

    colormatrix_halfspinor( (void*)C.D , (void*)Fmunu[i].O[6+2*mu] ,
			    (void*)S[ Sfwd2 ].D ) ;
    dagger_gauge( (void*)b , (void*)Fmunu[i].O[7+2*mu] ) ;
    colormatrix_halfspinor( (void*)D.D , (void*)b , (void*)S[ Sbck2 ].D ) ;

    // do the sum der = -5/2 S[i] + 4/3(A+B) - 1/12(C+D) 
    sq_inner( (void*)res.D , (void*)A.D , (void*)B.D ,
	      (void*)C.D , (void*)D.D , (void*)S[i].D ) ;
  }

  sigmaB_halfspinor( der , Fmunu[i] , res ) ;
  
  return ;
}
