/**
   @file mmul.c
   @brief matrix multiplies for colormatrices
 **/
#include "common.h"

#ifdef HAVE_EMMINTRIN_H

// simple NxN square matrix multiplication a = b.c
void 
multab( __m128d *__restrict a , 
	const __m128d *__restrict b , 
	const __m128d *__restrict c )
{
#if NC==3
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 8 ) ) ) ;
#elif NC == 2
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 3 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 3 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 2 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 3 ) , *( c + 3 ) ) ) ; 
#else
  int i , j , m ;
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m++  ) {
	sum = _mm_add_pd( sum , SSE2_MUL( *( b + m + NC*i ) , *( c + j + m*NC ) ) ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}

// 3x3 mult a = ( b^{\dagger} ).c 
void 
multabdag( __m128d *__restrict a , 
	   const __m128d *__restrict b , 
	   const __m128d *__restrict c )
{
#if NC == 3
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 0 ) ) ,
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 8 ) ) ) ;
#elif NC == 2
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 2 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 2 ) , *( c + 3 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 3 ) ) ) ;
#else
  int i , j , m ;
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = _mm_setzero_pd(  ) ;
      for( m = 0 ; m < NC ; m++ ) {
	sum = _mm_add_pd( sum , SSE2_MULCONJ( *( b + i + NC*m ) , *( c + j + NC*m ) ) ) ;
      }
      *a = sum ; a++ ;
    }
  }
#endif
  return ;
}

// a = b * c^{\dagger}
void 
multab_dag( __m128d *__restrict a , 
	    const __m128d *__restrict b , 
	    const __m128d *__restrict c )
{
#if NC == 3
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 6 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 7 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 8 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 6 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJ( *( b + 7 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 8 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 6 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJ( *( b + 7 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 8 ) , *( c + 8 ) ) ) ; a++ ;
#elif NC == 2
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 1 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 3 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 3 ) , *( c + 1 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 2 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJ( *( b + 3 ) , *( c + 3 ) ) ) ;
#else 
  int i , j , m ;
  register __m128d sum ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m++ ) {
	sum = _mm_add_pd( sum , SSE2_MUL_CONJ( *( b + m + NC*i ) , *( c + m + NC*j ) ) ) ;
      }
      *a = sum ; a++ ;
    }
  }
#endif
  return ;
}

// a = b^{\dagger} * c^{\dagger}
void 
multab_dagdag( __m128d *__restrict a , 
	       const __m128d *__restrict b , 
	       const __m128d *__restrict c )
{
#if NC == 3
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 8 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 8 ) ) ) ;
#elif NC == 2
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 1 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 3 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 1 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 3 ) ) ) ;
#else
  int i , j , m ;
  register __m128d sum ;
  const __m128d *pC , *pB ;
  for( i = 0 ; i < NC ; i++ ) {
    pC = c ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      sum = _mm_setzero_pd( ) ;
      for( m = 0 ; m < NC ; m++ ) {
	sum = _mm_add_pd( sum , SSE2_MUL_CONJCONJ( *( pB + i ) , *( pC + m ) ) ) ;
	pB += NC ;
      }
      *a = sum ; a++ ;
      pC += NC ;
    }
  }
#endif
  return ;
}

#endif
