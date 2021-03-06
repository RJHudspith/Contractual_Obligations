/**
   @file bar_ops_SSE.c
   @brief SSEd baryon operations
 */
#include "common.h"

#include "spinor_ops.h"  // spinor_zero_site()

#ifdef HAVE_EMMINTRIN_H

// This contracts the diquark with the remaining propagator
// This does the color trace Tr[ A . B^T ] 
double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const size_t d0 ,
		 const size_t d1 ,
		 const size_t d2 ,
		 const size_t d3 )
{
  const __m128d *d = (const __m128d*)DiQ.D[d0][d1].C ;
  const __m128d *s = (const __m128d*)S.D[d2][d3].C ;
#if NC == 3
  register __m128d sum ;
  sum = _mm_add_pd( SSE2_MUL( *d , *s ) , 
		    _mm_add_pd( SSE2_MUL( *(d+1) , *(s+1) ) , 
				SSE2_MUL( *(d+2) , *(s+2) ) ) ) ; 
  d+=3 ; s+=3 ;
  sum = _mm_add_pd( _mm_add_pd( sum , SSE2_MUL( *(d) , *(s) ) ) ,
		    _mm_add_pd( SSE2_MUL( *(d+1) , *(s+1) ) , 
				SSE2_MUL( *(d+2) , *(s+2) ) ) ) ; 
  d+=3 ; s+=3 ;
  sum = _mm_add_pd( _mm_add_pd( sum , SSE2_MUL( *(d) , *(s) ) ) ,
		    _mm_add_pd( SSE2_MUL( *(d+1) , *(s+1) ) , 
				SSE2_MUL( *(d+2) , *(s+2) ) ) ) ; 
#elif NC == 2
  register __m128d sum = SSE2_MUL( *d , *s) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
#else
  register __m128d sum = _mm_setzero_pd() ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  }
#endif
  double complex res ;
  _mm_store_pd( (void*)&res , sum ) ;
  return res ;
}

#if NC == 3
// cross color 
static inline void
cross_color( __m128d *__restrict a ,
	     const __m128d *__restrict b ,
	     const __m128d *__restrict c ,
	     const size_t id1 ,
	     const size_t id2 )
{
  // top
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 0 ] , c[ id2 + 0 ] ) ,
				    SSE2_MUL( b[ id2 + 0 ] , c[ id1 + 0 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 0 ] , c[ id2 + 3 ] ) ,
				    SSE2_MUL( b[ id2 + 0 ] , c[ id1 + 3 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 0 ] , c[ id2 + 6 ] ) ,
				    SSE2_MUL( b[ id2 + 0 ] , c[ id1 + 6 ] ) ) ) ;
  a++ ;
  // middle
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 3 ] , c[ id2 + 0 ] ) ,
				    SSE2_MUL( b[ id2 + 3 ] , c[ id1 + 0 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 3 ] , c[ id2 + 3 ] ) ,
				    SSE2_MUL( b[ id2 + 3 ] , c[ id1 + 3 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 3 ] , c[ id2 + 6 ] ) ,
				    SSE2_MUL( b[ id2 + 3 ] , c[ id1 + 6 ] ) ) ) ;
  a++ ;
  // bottom
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 6 ] , c[ id2 + 0 ] ) ,
				    SSE2_MUL( b[ id2 + 6 ] , c[ id1 + 0 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 6 ] , c[ id2 + 3 ] ) ,
				    SSE2_MUL( b[ id2 + 6 ] , c[ id1 + 3 ] ) ) ) ;
  a++ ;
  *a = _mm_add_pd( *a , _mm_sub_pd( SSE2_MUL( b[ id1 + 6 ] , c[ id2 + 6 ] ) ,
				    SSE2_MUL( b[ id2 + 6 ] , c[ id1 + 6 ] ) ) ) ;
  a++ ;
  return ;
}
#endif

// This carries out the color cross product and traces one set of Dirac indices.
// The result forms a diquark-type object
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S )
{
#if NC == 3
  // temporary 3-spinor space
  struct spinor T3SNK[ 3 ] ;
  // initialise to zero, apparently {} is a GNU extension
  spinor_zero_site( &T3SNK[0] ) ;
  spinor_zero_site( &T3SNK[1] ) ;
  spinor_zero_site( &T3SNK[2] ) ;
  // Sink cross color and trace, this leaves only one set of dirac indices
  size_t i , j , d ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      __m128d *a0 = (__m128d*)T3SNK[0].D[i][j].C ;
      __m128d *a1 = (__m128d*)T3SNK[1].D[i][j].C ;
      __m128d *a2 = (__m128d*)T3SNK[2].D[i][j].C ;
      // Dirac sink trace
      for( d = 0 ; d < NS ; d++ ) {
	const __m128d *b = (const __m128d*)S.D[i][d].C ;
	const __m128d *c = (const __m128d*)DiQ->D[j][d].C ;
	cross_color( a0 , b , c , 1 , 2 ) ;
	cross_color( a1 , b , c , 2 , 0 ) ;
	cross_color( a2 , b , c , 0 , 1 ) ; 
      }
    }
  }
  __m128d *D1 = (__m128d*)DiQ -> D ;
  const __m128d *t0 = (const __m128d*)T3SNK[0].D ;
  const __m128d *t1 = (const __m128d*)T3SNK[1].D ;
  const __m128d *t2 = (const __m128d*)T3SNK[2].D ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *D1 = _mm_sub_pd( *( t0 + NC + 2 ) , *( t0 + 1 + 2 * NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + NC + 2 ) , *( t1 + 1 + 2 * NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + NC + 2 ) , *( t2 + 1 + 2 * NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t0 + 2 * NC ) , *( t0 + 2 ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + 2 * NC ) , *( t1 + 2 ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + 2 * NC ) , *( t2 + 2 ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t0 + 1 ) , *( t0 + NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + 1 ) , *( t1 + NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + 1 ) , *( t2 + NC ) ) ; D1++ ;
    t0 += NCNC ; t1 += NCNC ; t2 += NCNC ;
  }
#elif NC == 2
  fprintf( stderr , "[CROSS COLOR TRACE] NC = %d not supported\n" , NC ) ;
  fprintf( stderr , "please call the diquark code instead of baryons\n" ) ;
  exit(1) ;
#else
  fprintf( "[CROSS COLOR TRACE] NC = %d not supported\n" , NC ) ;
  exit(1) ;
#endif
  return ;
}

#endif
