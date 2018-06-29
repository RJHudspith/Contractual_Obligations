/**
   @file grad_4.c
   @brief fourth order gradient terms
 */
#include "common.h"

#include "matrix_ops.h"     // add_mat()
#include "mmul.h"           // multab()
#include "halfspinor_ops.h" // zero_halfspinor() etc

// some little inline tricks here
#ifdef HAVE_IMMINTRIN_H

#if NC==3
#define inner_unroll_su3()\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\
  *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
  H++ ; A++ ; B++ ; E++ ; F++ ;						\

#define inner2_unroll_su3()\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\
  *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ;						\

#define grad4_inner_unroll()						\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\
  *H = _mm_add_pd( *H , _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,_mm_add_pd( _mm_add_pd( *C , *D ) , _mm_mul_pd( f2 , *S ) ) ) ) ; \
  H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;					\

#endif

static inline void
sqsq_inner( __m128d *H ,
	    const __m128d *A ,
	    const __m128d *B ,
	    const __m128d *E ,
	    const __m128d *F )
{
  register const __m128d m12 = _mm_set_pd( -12 , -12 ) ;
#if NC==3 && ND==4
  inner_unroll_su3() ;
  inner_unroll_su3() ;
  inner_unroll_su3() ;
  inner_unroll_su3() ;
#else
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H = _mm_add_pd( SSE2_FMA( m12 , _mm_add_pd( *A , *B ) , _mm_add_pd( *E , *F ) ) , *H ) ; \
    H++ ; A++ ; B++ ; E++ ; F++ ;	
  }
#endif
}

static inline void
sqsq_inner2( __m128d *H ,
	     const __m128d *A ,
	     const __m128d *B ,
	     const __m128d *C ,
	     const __m128d *D )
{
#if NC==3 && ND==4
  inner2_unroll_su3() ;
  inner2_unroll_su3() ;
  inner2_unroll_su3() ;
  inner2_unroll_su3() ;
#else
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H = _mm_add_pd( _mm_add_pd( _mm_add_pd( *A , *B ) , _mm_add_pd( *C , *D ) ) , * H ) ;
    H++ ; A++ ; B++ ; C++ ; D++ ; 	
  }
#endif
}

static inline void
grad4_inner( __m128d *H ,
	     const __m128d *A ,
	     const __m128d *B ,
	     const __m128d *C ,
	     const __m128d *D ,
	     const __m128d *S )
{
  register const __m128d f1 = _mm_setr_pd( -4 , -4 ) ;
  register const __m128d f2 = _mm_setr_pd(  6 ,  6 ) ;
#if NC==3 && ND==4
  grad4_inner_unroll() ;
  grad4_inner_unroll() ;
  grad4_inner_unroll() ;
  grad4_inner_unroll() ;
#else
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H = _mm_add_pd( *H ,
		     _mm_add_pd( _mm_mul_pd( f1 , _mm_add_pd( *A , *B ) ) ,
				 _mm_add_pd( _mm_add_pd( *C , *D ) , 
					     _mm_mul_pd( f2 , *S ) ) ) ) ;
    H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;
  }
#endif
  return ;
}

#else

static inline void
sqsq_inner( double complex *H ,
	    const double complex *A ,
	    const double complex *B ,
	    const double complex *E ,
	    const double complex *F )
{
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H += -12 * ( *A + *B ) + ( *E + *F ) ; 
    H++ ; A++ ; B++ ; E++ ; F++ ;
  }
  return ;
}

static inline void
sqsq_inner2( double complex *H ,
	     const double complex *A ,
	     const double complex *B ,
	     const double complex *C ,
	     const double complex *D )
{
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H += *A + *B + *C + *D ;
    H++ ; A++ ; B++ ; C++ ; D++ ; 	
  }
  return ;
}

static inline void
grad4_inner( double complex *H ,
	     const double complex *A ,
	     const double complex *B ,
	     const double complex *C ,
	     const double complex *D ,
	     const double complex *S )
{
  size_t i ;
  for( i = 0 ; i < ND*NCNC ; i++ ) {
    *H += -4 * (*A + *B ) + ( *C + *D ) + 6 * ( *S ) ;
    H++ ; A++ ; B++ ; C++ ; D++ ; S++ ;
  }
  return ;
}

#endif

// I make this
// U(x)U(x+\mu)S(x+\mu) - 4U(x)S(x+\mu) + 6S(x)
// - 4U^\dagger(x-mu)S(x-\mu) + U^\dagger(x-\mu)U^\dagger(x-2\mu)S(x-2\mu)
void
grad4( struct halfspinor *der ,
       const struct halfspinor *S ,
       const struct field *Fmunu ,
       const size_t i ,
       const size_t t ,
       const size_t mu )
{
#ifdef HAVE_IMMINTRIN_H
  __m128d C[ NCNC ] , D[ NCNC ] ;
#else
  double complex C[ NCNC ] , D[ NCNC ] ;
#endif
  
  const size_t Uidx = i + t*LCU ;
  
  const size_t Sfwd  = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
  
  const size_t Sbck  = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;
  
  const size_t Ubck = lat[ Uidx ].back[mu] ;

  zero_halfspinor( der ) ;
  dagger_gauge( C , (const void*)lat[ Ubck ].O[mu] ) ;
  dagger_gauge( D , (const void*)Fmunu[i].O[7+2*mu] ) ;

  struct halfspinor a , b , c , d ;
  colormatrix_halfspinor( (void*)a.D , (const void*)lat[ Uidx ].O[mu] ,
			  (const void*)S[Sfwd].D ) ;
  colormatrix_halfspinor( (void*)b.D , C , (const void*)S[Sbck].D ) ;
  colormatrix_halfspinor( (void*)c.D , (const void*)Fmunu[i].O[6+2*mu] ,
			  (const void*)S[Sfwd2].D ) ;
  colormatrix_halfspinor( (void*)d.D , D , (const void*)S[Sbck2].D ) ;

  grad4_inner( (void*)der -> D ,
	       (const void*)a.D ,
	       (const void*)b.D ,
	       (const void*)c.D ,
	       (const void*)d.D ,
	       (const void*)S[i].D ) ;
 
  return ;
}

// extra special \grad^2 \grad^2 code. This thing is a beast and pretty
// slow but it uses no temporaries
void
grad_sqsq( struct halfspinor *der2 ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t )
{
  struct halfspinor a , b , c , d ;
#ifdef HAVE_IMMINTRIN_H
  __m128d A[ NCNC ] , B[ NCNC ] , C[ NCNC ] , D[ NCNC ] ;
#else
  double complex A[ NCNC ] , B[ NCNC ] , C[ NCNC ] , D[ NCNC ] ;
#endif
  zero_halfspinor( der2 ) ;
  
  const size_t Uidx = i + t*LCU ;
  size_t mu , nu ;
  
  // beginning of the mu-nu loop
  for( mu = 0 ; mu < ND-1 ; mu++ ) {

    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    const size_t S_PmPm = lat[ lat[i].neighbor[mu] ].neighbor[mu] ;
    const size_t S_MmMm = lat[ lat[i].back[mu] ].back[mu] ;

    dagger_gauge( (void*)C , (const void*)lat[ Ubck ].O[mu] ) ;
    dagger_gauge( (void*)D , (const void*)Fmunu[i].O[7+2*mu]  ) ;

    colormatrix_halfspinor( (void*)a.D , (const void*)lat[ Uidx ].O[mu] , (const void*)S[ Sfwd ].D ) ;
    colormatrix_halfspinor( (void*)b.D , C , (const void*)S[ Sbck ].D ) ;
    colormatrix_halfspinor( (void*)c.D , (const void*)Fmunu[i].O[6+2*mu] , (const void*)S[ S_PmPm ].D ) ;
    colormatrix_halfspinor( (void*)d.D , D , (const void*)S[ S_MmMm ].D ) ;

    sqsq_inner( (void*)der2 -> D ,
		(const void*)a.D , (const void*)b.D ,
		(const void*)c.D , (const void*)d.D ) ;

    // inner nu loop
    for( nu = 0 ; nu < ND-1 ; nu++ ) {

      if( mu == nu ) continue ;
      
      const size_t S_PmPn = lat[ lat[i].neighbor[mu] ].neighbor[nu] ;
      const size_t S_PmMn = lat[ lat[i].neighbor[mu] ].back[nu] ;
      const size_t S_MmPn = lat[ lat[i].back[mu] ].neighbor[nu] ;
      const size_t S_MmMn = lat[ lat[i].back[mu] ].back[nu] ;

      const size_t U_Pm   = lat[ Uidx ].neighbor[mu] ;
      const size_t U_Mm   = lat[ Uidx ].back[mu] ;
      const size_t U_PmMn = lat[ lat[ Uidx ].neighbor[mu] ].back[nu] ;
      const size_t U_MmMn = lat[ lat[ Uidx ].back[mu] ].back[nu] ;
      
      // precompute these things
      multab( A , (void*)lat[ Uidx ].O[mu] , (void*)lat[ U_Pm ].O[nu] ) ;
      multabdag( B , (void*)lat[ U_Mm ].O[mu] , (void*)lat[ U_Mm ].O[nu] ) ;
      multab_dag( C , (void*)lat[ Uidx ].O[mu] , (void*)lat[ U_PmMn ].O[nu]  ) ;
      multab_dagdag( D , (void*)lat[ U_Mm ].O[mu] , (void*)lat[ U_MmMn ].O[nu] ) ;

      colormatrix_halfspinor( (void*)a.D , A , (const void*)S[ S_PmPn ].D ) ;
      colormatrix_halfspinor( (void*)b.D , B , (const void*)S[ S_MmPn ].D ) ;
      colormatrix_halfspinor( (void*)c.D , C , (const void*)S[ S_PmMn ].D ) ;
      colormatrix_halfspinor( (void*)d.D , D , (const void*)S[ S_MmMn ].D ) ;

      sqsq_inner2( (void*)der2 -> D ,
		   (const void*)a.D , (const void*)b.D ,
		   (const void*)c.D , (const void*)d.D ) ;
    }
  }
  
  // can add this right at the end for slightly fewer instruction requests
  // it is pretty funny that the answer to life the universe and everything
  // is in here
  halfspinor_Saxpy( der2 , S[i] , 42. ) ;

  return ;
}
