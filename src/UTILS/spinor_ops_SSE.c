/**
   @file spinor_ops.c
   @brief gauge*spinor, spinor*gauge and daggered variants

   SSE2 variants
 */
#include "common.h"

#include "matrix_ops.h" // NC*NC multiplies, daggers, sums and traces
#include "spinor_ops.h" // so that I can alphabetise

#ifdef HAVE_EMMINTRIN_H

// atomically add two spinors
static void
add_spinors( __m128d *SUM ,
	     const __m128d *S )
{
  int i ;
#if NC == 3 
  for( i = 0 ; i < NSNS ; i++ ) {
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
  int i ;
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
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *S = zero ; S++ ;
  }
#endif
}

// equate spinors
void
equate_spinor( void *S ,
	       const void *S2 )
{
  __m128d *s = (__m128d*)S ;
  const __m128d *s2 = (const __m128d*)S2 ;
  // probably be better calling out to memcpy ....
  int i ;
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
  int i ;
#if NC == 3 
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
    *s = -*s2 ; s++ ; s2++ ;
  }
#else
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s = -*s2 ; s++ ; s2++ ;
  }
#endif
  return ;
}

// flip a spinor, call with flipsign_spinor( (double complex*)S.D )
void
flipsign_spinor( void *S ) 
{
  int i ;
  __m128d *s = (__m128d*)S ;
#if NC == 3
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
    *s = -*s ; s++ ;
  }
#else
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
     *s = -*s ; s++ ;
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
  int d1d2 ;
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
  int d1d2 ;
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
  int d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab_dag( r , l , s ) ;
    r += NCNC ; s += NCNC ;
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
  int d1d2 ;
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
  int d1d2 ;
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
  int d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    multab_dag( r , s , l ) ;
    r += NCNC ; s += NCNC ;
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
  int i ;
  for( i = 0 ; i < VOL3 ; i++ ) {
    add_spinors( tSUM , tS ) ; tS += NSNS*NCNC ;
  }
  return ;
}

// multiplies two spinors A = B * A
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B )
{
  struct spinor tmp = *A ;
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      // color matrix multiply
      multab( (__m128d*)A -> D[ d1 ][ d2 ].C , 
	      (const __m128d*)B.D[ d1 ][ d2 ].C , 
	      (const __m128d*)tmp.D[ d2 ][ d1 ].C ) ;
    }
  }
  return ;
}

// multithreaded zero a spinor over a timeslice
void
spinor_zero( void *S )
{
  __m128d *s = (__m128d*)S ;
  int i ;
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

#endif
