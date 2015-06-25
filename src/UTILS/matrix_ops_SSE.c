/**
   @file matrix_ops_SSE.c
   @brief simple matrix multiplies (SSEd version)
 */

#include "common.h"
#include "matrix_ops.h"

#ifdef HAVE_EMMINTRIN_H

// add two color matrices
void
add_mat( __m128d *__restrict a ,
	 const __m128d *b )
{
#if NC == 3
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
#elif NC == 2
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *a += *b ; a++ ; b++ ;
  }
#endif
}

// equate two color matrices
// equate two color matrices
void
colormatrix_equiv( double complex *__restrict a ,
		   const double complex *__restrict b )
{
#if NC == 3
  a[0] = b[0] ; a[1] = b[1] ; a[2] = b[2] ; 
  a[3] = b[3] ; a[4] = b[4] ; a[5] = b[5] ; 
  a[6] = b[6] ; a[7] = b[7] ; a[8] = b[8] ; 
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = b[i] ;
  }
#endif
  return ;
}

// float matrix to double
void
colormatrix_equiv_f2d( double complex a[ NCNC ] ,
		       const float complex b[ NCNC ] )
{
#if NC == 3
  a[0] = (double complex)b[0] ; a[1] = (double complex)b[1] ; a[2] = (double complex)b[2] ; 
  a[3] = (double complex)b[3] ; a[4] = (double complex)b[4] ; a[5] = (double complex)b[5] ; 
  a[6] = (double complex)b[6] ; a[7] = (double complex)b[7] ; a[8] = (double complex)b[8] ; 
#elif NC == 2
  a[0] = (double complex)b[0] ; a[1] = (double complex)b[1] ;
  a[2] = (double complex)b[2] ; a[3] = (double complex)b[3] ;
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = (double complex)b[i] ;
  }
#endif
  return ;
}

// trace of a color matrix
__m128d
colortrace( const __m128d *a )
{
#if NC == 3
  return _mm_add_pd( _mm_add_pd( a[ 0 ] , a[ 4 ] ) , a[8] ) ; 
#else
  register __m128d sum = _mm_setzero_pd( ) ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    sum = _mm_add_pd( sum , a[ i*( NC + 1 ) ] ) ;
  }
  return sum ;
#endif
}

// is just Tr( a * b )
__m128d
colortrace_prod( const __m128d *a ,
		 const __m128d *b )
{
#if NC == 3
  register __m128d sum = _mm_add_pd( SSE2_MUL( a[0] , b[0] ) , SSE2_MUL( a[1] , b[3] ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( SSE2_MUL( a[2] , b[6] ) , SSE2_MUL( a[3] , b[1] ) ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( SSE2_MUL( a[4] , b[4] ) , SSE2_MUL( a[5] , b[7] ) ) ) ;
  sum = _mm_add_pd( sum , _mm_add_pd( SSE2_MUL( a[6] , b[2] ) , SSE2_MUL( a[7] , b[5] ) ) ) ;
  return _mm_add_pd( sum , SSE2_MUL( a[8] , b[8] ) ) ;
#elif NC == 2
  register __m128d sum = _mm_add_pd( SSE2_MUL( a[0] , b[0] ) , SSE2_MUL( a[1] , b[2] ) ) ;
  return _mm_add_pd( sum , _mm_add_pd( SSE2_MUL( a[2] , b[1] ) , SSE2_MUL( a[3] , b[3] ) ) ) ;
#else
  register __m128d sum = _mm_setzero_pd( ) ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    int j ;
    for( j = 0 ; j < NC ; j++ ) {
      sum = _mm_add_pd( sum , SSE2_MUL( a[j+i*NC] , b[i+j*NC] ) ) ;
    }
  }
  return sum ;
#endif
}

// does res = constant * U
void
constant_mul_gauge( double complex *__restrict res , 
		    const double complex constant ,
		    const double complex *__restrict U ) 
{
#if NC == 3
  res[0] = constant * U[0] ; res[1] = constant * U[1] ; res[2] = constant * U[2] ;
  res[3] = constant * U[3] ; res[4] = constant * U[4] ; res[5] = constant * U[5] ;
  res[6] = constant * U[6] ; res[7] = constant * U[7] ; res[8] = constant * U[8] ;
#elif NC == 2
  res[0] = constant * U[0] ; res[1] = constant * U[1] ; 
  res[2] = constant * U[2] ; res[3] = constant * U[3] ;
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    res[ i ] = constant * U[ i ] ;
  }
#endif
  return ;
}

// daggers the matrix U into res
void
dagger_gauge( __m128d *__restrict res ,
	      const __m128d *__restrict U )
{
#if NC == 3
  *res = SSE2_CONJ( *( U + 0 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 3 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 6 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 1 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 4 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 7 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 2 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 5 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 8 ) ) ; res++ ; 
#elif NC == 2
  *res = SSE2_CONJ( *( U + 0 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 2 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 1 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 3 ) ) ; res++ ; 
#else
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      *res = SSE2_CONJ( U[ i + j * NC ] ) ; res++ ;
    }
  }
#endif
  return ;
}

//  the column-pivoted LU decomposition determinant
//  does not save L, just need the diagonal of U as determinant is product
//  of these elements
double complex
LU_det( const int N , 
	const double complex U[ N*N ] )
{
  int i , j , l , piv , perms = 0 ;
  double complex a[ N*N ] , dt , determinant = 1. ;
  // workspace is double precision
  for( i = 0 ; i < N*N ; i++ ) {
    a[ i ] = (double complex)U[ i ] ;
  }
  double attempt , best ; 
  for( i = 0 ; i < N-1 ; i++ ) {
    // swap rows s.t a[i] has largest pivot number first
    best = creal( a[i*(N+1)] ) * creal( a[i*(N+1)] ) 
         + cimag( a[i*(N+1)] ) * cimag( a[i*(N+1)] ) ;
    piv = i ;
    // again only care about the pivots below i
    for( j = i+1 ; j < N ; j++ ) {
      attempt = creal( a[i+j*N] ) * creal( a[i+j*N] ) 
	      + cimag( a[i+j*N] ) * cimag( a[i+j*N] ) ;
      if( attempt > best ) { 
	piv = j ; 
	best = attempt ; 
      }
    }
    if( a[i+piv*N] == 0.0 ) { 
      printf( "[DETERMINANT] LU  Singular Matrix!!!\n" ) ;
      return 0.0 ;
    }
    if( piv != i ) {
      // swap rows
      for( l = 0 ; l < N ; l++ ) {
	dt         = a[l+i*N] ;
	a[l+i*N]   = a[l+piv*N] ;
	a[l+piv*N] = dt ;
      }
      perms++ ;
    }
    // perform gaussian elimination
    dt = 1.0 / a[ i*(N+1) ] ;
    double complex *pA = a + i*N ;
    for( j = N-1 ; j > i ; j-- ) { // go up in other column
      register double complex fac1 = a[ i + j*N ] * dt ; 
      // go along the row performing the subtraction, there is no point in
      // subtracting elements where we have determined the best pivot, just the
      // columns to the right of the pivot
      for( l = i + 1 ; l < N ; l++ ) {
	a[ l + j*N ] -= creal( fac1 ) * creal( pA[l] ) - cimag( fac1 ) * cimag( pA[l] ) 
	        + I * ( creal( fac1 ) * cimag( pA[l] ) + cimag( fac1 ) * creal( pA[l] ) ) ;
      }
    }
    determinant *= a[ i*(N+1) ] ;
  }
  determinant *= a[ N*N-1 ] ;
  return perms&1 ? -determinant : determinant ;
}

// simple NxN square matrix multiplication a = b.c
void 
multab( __m128d *__restrict a , 
	const __m128d *__restrict b , 
	const __m128d *__restrict c )
{
#if NC==3
  // first row
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 2 ) , *( c + 8 ) ) ) ; a++ ;
  // second row
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 3 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 4 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 5 ) , *( c + 8 ) ) ) ; a++ ;
  // third row
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL( *( b + 6 ) , *( c + 2 ) ) , 
		   SSE2_MUL( *( b + 7 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL( *( b + 8 ) , *( c + 8 ) ) ) ; a++ ;
#elif NC == 2
  // a[0] = b[0] * c[0] + b[1] * c[2]
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 2 ) ) ) ; a++ ;
  // a[1] = b[0] * c[1] + b[1] * c[3]
  *a = _mm_add_pd( SSE2_MUL( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 1 ) , *( c + 3 ) ) ) ; a++ ;
  // a[2] = b[2] * c[0] + b[3] * c[2]
  *a = _mm_add_pd( SSE2_MUL( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL( *( b + 3 ) , *( c + 2 ) ) ) ; a++ ;
  // a[3] = b[2] * c[1] + b[3] * c[3]
  *a = _mm_add_pd( SSE2_MUL( *( b + 2 ) , *( c + 1 ) ) , 
		   SSE2_MUL( *( b + 3 ) , *( c + 3 ) ) ) ; a++ ;
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
  // first row
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 6 ) , *( c + 8 ) ) ) ; a++ ;
  // second
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 4 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 7 ) , *( c + 8 ) ) ) ; a++ ;
  // third
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 3 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 6 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 7 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 2 ) , *( c + 2 ) ) , 
		   SSE2_MULCONJ( *( b + 5 ) , *( c + 5 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MULCONJ( *( b + 8 ) , *( c + 8 ) ) ) ; a++ ;
#elif NC == 2
  // a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2]
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 2 ) , *( c + 2 ) ) ) ; a++ ;
  // a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3]
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 0 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 2 ) , *( c + 3 ) ) ) ; a++ ;
  // a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2]
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 2 ) ) ) ; a++ ;
  // a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3]
  *a = _mm_add_pd( SSE2_MULCONJ( *( b + 1 ) , *( c + 1 ) ) , 
		   SSE2_MULCONJ( *( b + 3 ) , *( c + 3 ) ) ) ; a++ ;
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
  // first row
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 2 ) , *( c + 8 ) ) ) ; a++ ;
  // second
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 3 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJ( *( b + 4 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJ( *( b + 5 ) , *( c + 8 ) ) ) ; a++ ;
  // third
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
  //a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] )
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 1 ) ) ) ; a++ ;
  //a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) 
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJ( *( b + 1 ) , *( c + 3 ) ) ) ; a++ ;
  //a[2] = b[2] * conj( c[0] ) + b[3] * conj( c[1] ) 
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJ( *( b + 3 ) , *( c + 1 ) ) ) ; a++ ;
  //a[3] = b[2] * conj( c[2] ) + b[3] * conj( c[3] )
  *a = _mm_add_pd( SSE2_MUL_CONJ( *( b + 2 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJ( *( b + 3 ) , *( c + 3 ) ) ) ; a++ ;
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
  // first row
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 6 ) , *( c + 8 ) ) ) ; a++ ;
  // second
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 4 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 7 ) , *( c + 8 ) ) ) ; a++ ;
  // third
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 1 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 2 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 3 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 4 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 5 ) ) ) ; a++ ;
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 6 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 5 ) , *( c + 7 ) ) ) ;
  *a = _mm_add_pd( *a , SSE2_MUL_CONJCONJ( *( b + 8 ) , *( c + 8 ) ) ) ; a++ ;
#elif NC == 2
  //a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] )
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 1 ) ) ) ; a++ ;
  //a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] )
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 0 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 2 ) , *( c + 3 ) ) ) ; a++ ;
  //a[2] = conj( b[1] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] )
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 0 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 1 ) ) ) ; a++ ;
  //a[3] = conj( b[1] ) * conj( c[2] ) + conj( b[3] ) * conj( c[3] )
  *a = _mm_add_pd( SSE2_MUL_CONJCONJ( *( b + 1 ) , *( c + 2 ) ) , 
		   SSE2_MUL_CONJCONJ( *( b + 3 ) , *( c + 3 ) ) ) ; a++ ;
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

// print "a" to stdout
void
print_colormatrix( const double complex a[ NCNC ] )
{
  int i , j ;
  printf( "\n" ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      printf( "%1.3f %1.3f " , creal( a[j+i*NC] ) , cimag( a[j+i*NC] ) ) ;
    }
    printf( "\n" ) ;
  }
  return ;
}

#endif
