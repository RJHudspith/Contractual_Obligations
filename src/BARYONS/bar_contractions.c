/**
   @file bar_contractions.c
   @brief baryon contraction codes
 */

#include "common.h"

#include "spinor_ops.h"   // spinor_zero_site()

// This contracts the diquark with the remaining propagator
// This does the color trace Tr[ A . B ] 
const double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const int d0 ,
		 const int d1 ,
		 const int d2 ,
		 const int d3 )
{
#if (defined HAVE_IMMINTRIN_H) && (NC == 3)
  const __m128d *d = (const __m128d*)DiQ.D[d1][d0].C ;
  const __m128d *s = (const __m128d*)S.D[d2][d3].C ;
  register __m128d sum = _mm_setzero_pd() ;
  sum = SSE2_MUL( *d , *s) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  double complex res ;
  _mm_store_pd( (void*)&res , sum ) ;
  return res ;
#else
  int c1, c2 ;
  register double corrr = 0.0 , corri = 0.0 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      corrr += creal( DiQ.D[d1][d0].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] ) 
	- cimag( DiQ.D[d1][d0].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] );
      corri += creal( DiQ.D[d1][d0].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] ) 
	+ cimag( DiQ.D[d1][d0].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] );
    }
  }
  return corrr + I * corri;
#endif
}

// This carries out the color cross product and traces one set of Dirac indices.
// The results forms a diquark-type object
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S )
{
  // temporary 3-spinor space
  struct spinor T3SNK[ 3 ] ;
  spinor_zero_site( &T3SNK[ 0 ] ) ;
  spinor_zero_site( &T3SNK[ 1 ] ) ;
  spinor_zero_site( &T3SNK[ 2 ] ) ;

  int c1, c2, i, j, d;
  // Sink cross color and trace, this leaves only one set of dirac indices
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      // Dirac sink trace
      for( d = 0 ; d < NS ; d++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    T3SNK[0].D[i][j].C[c1][c2] += S.D[i][d].C[c1][1] * DiQ -> D[j][d].C[c2][2] ; 
	    T3SNK[0].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][2] * DiQ -> D[j][d].C[c2][1] ;
	    T3SNK[1].D[i][j].C[c1][c2] += S.D[i][d].C[c1][2] * DiQ -> D[j][d].C[c2][0] ;
	    T3SNK[1].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][0] * DiQ -> D[j][d].C[c2][2] ;
	    T3SNK[2].D[i][j].C[c1][c2] += S.D[i][d].C[c1][0] * DiQ -> D[j][d].C[c2][1] ;
	    T3SNK[2].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][1] * DiQ -> D[j][d].C[c2][0] ;
	  }
	}
      }
    }
  }

  int c3 ;
  // Source cross color
  for( c3 = 0 ; c3 < 3 ; c3++ ) {
    for( i = 0 ; i < NS ; i++ ) {
      for( j = 0 ; j < NS ; j++ ) {
	DiQ -> D[i][j].C[0][c3] = T3SNK[c3].D[i][j].C[1][2] - T3SNK[c3].D[i][j].C[2][1] ; 
	DiQ -> D[i][j].C[1][c3] = T3SNK[c3].D[i][j].C[2][0] - T3SNK[c3].D[i][j].C[0][2] ; 
	DiQ -> D[i][j].C[2][c3] = T3SNK[c3].D[i][j].C[0][1] - T3SNK[c3].D[i][j].C[1][0] ; 
      }            
    }
  }
  return ;
}
