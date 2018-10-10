/**
   @file spinor_ops_SSE.c
   @brief general spinor operations (SSEd version)
 */
#include "common.h"

#include "matrix_ops.h" // daggers, sums and traces
#include "mmul.h"       // NC*NC multiplies
#include "spinor_ops.h" // so that I can alphabetise

#ifdef HAVE_EMMINTRIN_H

// sum a propagator
static inline void
sum_spinor( __m128d *SUM ,
	    const __m128d *S )
{
  size_t i ;
#if ND == 4 
  for( i = 0 ; i < NCNC ; i++ ) {
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
  }
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *SUM = _mm_add_pd( *SUM , *S ) ; SUM++ ; S++ ;
  }
#endif
  return ;
}

// zero a spinor
static void
zero_spinor( __m128d *S )
{
  register const __m128d zero = _mm_setzero_pd( ) ;
  size_t i ;
#if NC == 3 
  for( i = 0 ; i < NSNS ; i++ ) {
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
  }
#elif NC == 2
  for( i = 0 ; i < NSNS ; i++ ) {
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
    *S = zero ; S++ ;
  }
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *S = zero ; S++ ;
  }
#endif
}

// atomically add spinors
void
add_spinors( struct spinor *A ,
	     const struct spinor B )
{
  sum_spinor( (__m128d*)A -> D , (const __m128d*)B.D ) ;
  return ;
}

// colortrace our spinor into a dirac matrix
void
colortrace_spinor( void *S1 ,
		   const void *S2 )
{
  __m128d *s = (__m128d*)S1 ;
  const __m128d *s2 = (const __m128d*)S2 ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = _mm_setzero_pd( ) ;
    #if NC == 3
    *s = _mm_add_pd( *s , *( s2 + 0 ) ) ; 
    *s = _mm_add_pd( *s , *( s2 + 4 ) ) ; 
    *s = _mm_add_pd( *s , *( s2 + 8 ) ) ; 
    #elif NC == 2
    *s = _mm_add_pd( *s , *( s2 + 0 ) ) ; 
    *s = _mm_add_pd( *s , *( s2 + 3 ) ) ; 
    #else
    size_t c ;
    for( c = 0 ; c < NC ; c++ ) {
      *s = _mm_add_pd( *s , *( s2 + c * ( NC + 1 ) ) ) ; 
    }
    #endif
    s2 += NCNC ;
    s++ ;
  }
  return ;
}

// equate spinors
void
equate_spinor( void *S ,
	       const void *S2 )
{
  __m128d *s = (__m128d*)S ;
  const __m128d *s2 = (const __m128d*)S2 ;
  // probably be better calling out to memcpy ....
  size_t i ;
#if NC == 3 
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
  }
#elif NC == 2
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
    *s = *s2 ; s++ ; s2++ ;
  }
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s = *s2 ; s++ ; s2++ ;
  }
#endif
  return ;
}

// equate one spinor to the minus of another
void
equate_spinor_minus( void *mS ,
		     const void *S )
{
  __m128d *s = (__m128d*)mS ;
  const __m128d *s2 = (const __m128d*)S ;
  size_t i ;
#if NC == 3 
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
  }
#elif NC == 2
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
  }
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s = SSE_FLIP( *s2 ) ; s++ ; s2++ ;
  }
#endif
  return ;
}

// flip a spinor, call with flipsign_spinor( (double complex*)S.D )
void
flipsign_spinor( void *S ) 
{
  size_t i ;
  __m128d *s = (__m128d*)S ;
#if NC == 3
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
  }
#elif NC == 2
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
    *s = SSE_FLIP( *s ) ; s++ ;
  }
#else
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    *s = SSE_FLIP( *s ) ; s++ ;
  }
#endif
  return ;
}

// multiply by a link :: res = ( link * S )
void
gauge_spinor( struct spinor *__restrict res ,
	      const double complex *__restrict link ,
	      const struct spinor S )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab( r , l , s ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// multiply by a daggered link res = link^{\dagger} * S
void
gaugedag_spinor( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multabdag( r , l , s ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// right multiply link by a daggered spinor res = link * S^{\dagger}
void
gauge_spinordag( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab_dag( r , l , s ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// set a spinor to the identity
void
identity_spinor( struct spinor *__restrict res )
{
  __m128d *r = (__m128d*)res -> D ;
  zero_spinor( r ) ;
  const register __m128d one = _mm_setr_pd( 1.0 , 0.0 ) ;
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    #if NC == 3
    *r = one ; *( r + NC + 1 ) = one ; *( r + 2 * ( NC + 1 ) ) = one ;
    #else
    size_t c ;
    for( c = 0 ; c < NC ; c++ ) {
      *( r + c * ( NC + 1 ) ) = one ;
    }
    #endif
    r += NCNC * ( NS + 1 ) ;
  }
  return ;
}

// multiply by a link :: res = S * link
void
spinor_gauge( struct spinor *__restrict res ,  
	      const struct spinor S ,
	      const double complex link[ NCNC ] )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab( r , s , l ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// multiply by a daggered link res = S^{\dagger} link
void
spinordag_gauge( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multabdag( r , s , l ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// right multiply by a daggered link res = S * link^{\dagger}
void
spinor_gaugedag( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] )
{
  __m128d *r = (__m128d*)res -> D ;
  const __m128d *l = (const __m128d*)link ;
  const __m128d *s = (const __m128d*)S.D ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab_dag( r , s , l ) ;
    r += NCNC ; s += NCNC ;
  }
  return ;
}

// multithreaded zero a spinor over a timeslice
void
spinor_zero( void *S )
{
  __m128d *s = (__m128d*)S ;
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    zero_spinor( s + i * ( NSNS * NCNC ) ) ;
  }
  return ;
}

// zero a spinor at a site
void
spinor_zero_site( void *S )
{
  zero_spinor( (__m128d*)S ) ;
}

// multiplies two spinors A = B * A
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B )
{
  struct spinor tmp = *A ;
  size_t d1 , d2 , d3 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      __m128d link[ NCNC ] , sum[ NCNC ] ;
      // zero
      for( d3 = 0 ; d3 < NCNC ; d3++ ) { 
	sum[ d3 ] = _mm_setzero_pd( ) ;
      }
      // color matrix multiply
      for( d3 = 0 ; d3 < NS ; d3++ ) { 
	multab( link , 
		(const __m128d*)B.D[ d1 ][ d3 ].C , 
		(const __m128d*)tmp.D[ d3 ][ d2 ].C ) ;
	add_mat( sum , link ) ;
      }
      colormatrix_equiv( (double complex*)A -> D[d1][d2].C ,
			 (const double complex*)sum ) ;
    }
  }
  return ;
}

// multiplies two spinors A = A * B
void
spinmul_atomic_right( struct spinor *A ,
		      const struct spinor B )
{
  struct spinor tmp = *A ;
  size_t d1 , d2 , d3 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      __m128d link[ NCNC ] , sum[ NCNC ] ;
      // zero
      for( d3 = 0 ; d3 < NCNC ; d3++ ) { 
	sum[ d3 ] = _mm_setzero_pd( ) ;
      }
      // color matrix multiply
      for( d3 = 0 ; d3 < NS ; d3++ ) { 
	multab( link , 
		(const __m128d*)tmp.D[ d1 ][ d3 ].C , 
		(const __m128d*)B.D[ d3 ][ d2 ].C ) ;
	add_mat( sum , link ) ;
      }
      colormatrix_equiv( (double complex*)A -> D[d1][d2].C ,
			 (const double complex*)sum ) ;
    }
  }
  return ;
}

// atomically add spinors
void
spinor_Saxpy( struct spinor *A ,
	      const double S ,
	      const struct spinor B )
{
  __m128d *s1 = (__m128d*)A -> D ;
  const __m128d *s2 = (const __m128d*)B.D ;
  register const __m128d s = _mm_setr_pd( S , S ) ;
  size_t i ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s1 = SSE2_FMA( *s2 , s , *s1 ) ; s1++ ; s2++ ;
  }
  return ;
}


// trace out our dirac indices
void
spintrace( void *S ,
	   const void *S2 )
{
  __m128d *s = (__m128d*)S ;
  const __m128d *s2 = (const __m128d*)S2 ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    #if NS == 4
    *s = *( s2 ) ;
    *s = _mm_add_pd( *s , *( s2 + NCNC * ( NS + 1 ) ) ) ;
    *s = _mm_add_pd( *s , *( s2 + 2 * NCNC * ( NS + 1 ) ) ) ;
    *s = _mm_add_pd( *s , *( s2 + 3 * NCNC * ( NS + 1 ) ) ) ;
    #else
    size_t j ;
    *s = _mm_setzero_pd( ) ;
    for( j = 0 ; j < NS ; j++ ) {
      *s = _mm_add_pd( *s , *( s2 + j * NCNC * ( NS + 1 ) ) ) ;
    }
    #endif
    s++ ;
    s2 ++ ;
  }
}

// atomically add spinors
void
sub_spinors( struct spinor *A ,
	     const struct spinor B )
{
  __m128d *s1 = (__m128d*)A -> D ;
  const __m128d *s2 = (const __m128d*)B.D ;
  size_t i ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s1 = _mm_sub_pd( *s1 , *s2 ) ; s1++ ; s2++ ;
  }
  return ;
}

// sums a propagator over a timeslice
void
sumprop( void *SUM ,
	 const void *S )
{
  __m128d *tSUM = (__m128d*)SUM ;
  zero_spinor( tSUM ) ;
  const __m128d *tS = (const __m128d*)S ;
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    sum_spinor( tSUM , tS ) ; tS += NSNS*NCNC ;
  }
  return ;
}

// sums propagators over spatial volume into spinor "SUM"
void
sumwalls( struct spinor *SUM ,
	  const struct spinor **S ,
	  const size_t Nprops )
{
  size_t mu ;
#pragma omp parallel for private(mu)
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    sumprop( &SUM[mu] , S[mu] ) ;
  }
  return ;
}

// dirac index transpose a spinor, returns S^T on stack
struct spinor
transpose_spinor( const struct spinor S )
{
  struct spinor ST ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    const size_t d1 = d1d2 / NS ;
    const size_t d2 = d1d2 % NS ;
    colormatrix_equiv( (void*)ST.D[d2][d1].C ,
		       (const void*)S.D[d1][d2].C ) ;
  }
  return ST ;
}

#endif
