/**
   @file grad.c
   @brief single derivatives
 */
#include "common.h"

#include "halfspinor_ops.h"  // Fmunu halfspinor
#include "matrix_ops.h"      // add_mat
#include "mmul.h"            // multab

#ifdef HAVE_IMMINTRIN_H

#if NC==3
#define unroll_improved_fac()\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;   \
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ; \
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;   \
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\
  *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;\
  A++ ; B++ ; C++ ; D++ ;\

#endif

static inline void
improved_fac( __m128d *A ,
	      const __m128d *B ,
	      const __m128d *C ,
	      const __m128d *D ,
	      const double U_0 )
{
  register const __m128d f2 = _mm_set_pd( 1/12. , 1/12. ) ;
  #ifdef LEGACY_NRQCD_COMPARE
  register const __m128d f1 = _mm_set_pd( 7.+1/(U_0*U_0) , 7.+1/(U_0*U_0) ) ;
  #else
  register const __m128d f1 = _mm_set_pd( 8. , 8. ) ;
  #endif
#if NC==3
  unroll_improved_fac() ;
  unroll_improved_fac() ;
  unroll_improved_fac() ;
  unroll_improved_fac() ;
  #else
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *A = _mm_mul_pd( f2 , SSE2_FMA( f1 , _mm_sub_pd( *A , *B ) , _mm_sub_pd( *D , *C ) ) ) ;
    A++ ; B++ ; C++ ; D++ ;
  }
  #endif
}

#else

static inline void
improved_fac( double complex *A ,
	      const double complex *B ,
	      const double complex *C ,
	      const double complex *D ,
	      const double U_0 )
{
  register const double f2 = 1/12.  ;
  #ifdef LEGACY_NRQCD_COMPARE
  register const double f1 = 7.+1/(U_0*U_0) ;
  #else
  register const double f1 = 8. ;
  #endif
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *A = f2 * ( f1 * ( *A - *B ) + ( *D - *C ) ) ;
    A++ ; B++ ; C++ ; D++ ;
  }
}

#endif

// computes improved derivative of S and puts it in der
void
grad_imp_LCU( struct halfspinor *der ,
	      const struct halfspinor *S ,
	      const struct field *Fmunu ,
	      const double U_0 ,
	      const size_t t ,
	      const size_t mu )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // some temporaries we need
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    
    const size_t Uidx = i + t*LCU ;    
    const size_t Ubck = lat[ Uidx ].back[mu] ;
    
    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      // computes der = U(x) S(x+\mu)
      multab( (void*)der[ i ].D[d] ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      // computes A = U^\dagger (x-\mu ) S(x-\mu)
      multabdag( (void*)A ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      // computes der = -( der - A )/2
      #ifdef LEGACY_NRQCD_COMPARE
      colormatrix_Sa_xmy( der[i].D[d] , A , (7+1/(U_0*U_0))/12. ) ;
      #else
      colormatrix_Sa_xmy( der[i].D[d] , A , 2./3. ) ;
      #endif
      
      // compute U(x)U(x+\mu)S(x+2\mu)
      multab( (void*)A ,
	      (void*)Fmunu[i].O[6+2*mu] ,
	      (void*)S[ Sfwd2 ].D[d] ) ;
      // compute U(x-mu)U(x-2mu)S(x-2\mu)
      multabdag( (void*)B ,
		 (void*)Fmunu[i].O[7+2*mu],
		 (void*)S[ Sbck2 ].D[d] ) ;
      colormatrix_Sa_xmy( A , B , -1./12. ) ;

      add_mat( (void*)der[i].D[d] , (void*)A ) ;
    }
  }
  return ;
}

// computes improved derivative of Fmunu.S and puts it in der
void
grad_imp_FMUNU( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  struct halfspinor A , B , C ;
#ifdef HAVE_IMMINTRIN_H
  __m128d D[ NCNC ] , E[ NCNC ] , F[ NCNC ] , G[ NCNC ] ;
#else
  double complex D[ NCNC ] , E[ NCNC ] , F[ NCNC ] , G[ NCNC ] ;
#endif
  
  const size_t Uidx  = i + t*LCU ;    
  const size_t Ubck  = lat[ Uidx ].back[mu] ;
  const size_t Sfwd  = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
  const size_t Sbck  = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;

  // precompute some common terms here -> products of links and clovers
  multab( D , (void*)lat[ Uidx ].O[mu] ,
	  (void*)Fmunu[ Sfwd ].O[ Fmunu_idx ] ) ;
  multabdag( E , (void*)lat[ Ubck ].O[mu] ,
	     (void*)Fmunu[ Sbck ].O[ Fmunu_idx ] ) ;
  multab( F , (void*)Fmunu[i].O[6+2*mu] ,
	  (void*)Fmunu[ Sfwd2 ].O[ Fmunu_idx ] ) ;
  multabdag( G , (void*)Fmunu[i].O[7+2*mu] ,
	     (void*)Fmunu[ Sbck2 ].O[ Fmunu_idx ] ) ;

  colormatrix_halfspinor( (void*)der -> D , D , (const void*)S[ Sfwd ].D ) ;
  colormatrix_halfspinor( (void*)A.D   , E , (const void*)S[ Sbck ].D ) ;
  colormatrix_halfspinor( (void*)B.D   , F , (const void*)S[ Sfwd2 ].D ) ;
  colormatrix_halfspinor( (void*)C.D   , G , (const void*)S[ Sbck2 ].D ) ;

  improved_fac( (void*)der -> D , (const void*)A.D ,
		(const void*)B.D , (const void*)C.D , U_0 ) ;

  return ;
}

// computes the improved derivative of S and multiplies on the left by Fmunu
void
FMUNU_grad_imp( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  // some temporaries we need
  struct halfspinor A , B , E , res ;
#ifdef HAVE_IMMINTRIN_H
  __m128d C[ NCNC ] , D[ NCNC ] ; 
#else
  double complex C[ NCNC ] , D[ NCNC ] ;
#endif
  const size_t Uidx = i + t*LCU ;    
  const size_t Ubck = lat[ Uidx ].back[mu] ;
    
  const size_t Sfwd = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
  const size_t Sbck = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;

  // it turns out a little faster to do this than call multabdag
  dagger_gauge( C , (void*)lat[ Ubck ].O[mu] ) ;
  dagger_gauge( D , (void*)Fmunu[i].O[7+2*mu] ) ;

  colormatrix_halfspinor( (void*)res.D ,
			  (const void*)lat[ Uidx ].O[mu] ,
			  (const void*)S[ Sfwd ].D ) ;
  colormatrix_halfspinor( (void*)A.D   , C , (const void*)S[ Sbck ].D ) ;
  colormatrix_halfspinor( (void*)B.D   ,
			  (const void*)Fmunu[i].O[6+2*mu] ,
			  (const void*)S[ Sfwd2 ].D ) ;
  colormatrix_halfspinor( (void*)E.D   , D , (const void*)S[ Sbck2 ].D ) ;

  improved_fac( (void*)res.D , (const void*)A.D ,
		(const void*)B.D , (const void*)E.D , U_0 ) ;

  // res is the improved gradient and we left multiply by the gauge field
  colormatrix_halfspinor( (void*)der -> D   ,
			  (const void*)Fmunu[i].O[ Fmunu_idx ] ,
			  (const void*)res.D ) ;
  return ;
}
