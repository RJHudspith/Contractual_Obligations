/**
   @file matrix_ops.c
   @brief simple matrix multiplies (SSE2 version)
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
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = (double complex)b[i] ;
  }
#endif
  return ;
}

// is just Tr( a * b )
__m128d
colortrace_prod( const __m128d *a ,
		 const __m128d *b )
{
  return 
    _mm_add_pd( SSE2_MUL( a[0] , b[0] ) , SSE2_MUL( a[1] , b[3] ) ) +
    _mm_add_pd( SSE2_MUL( a[2] , b[6] ) , SSE2_MUL( a[3] , b[1] ) ) +
    _mm_add_pd( SSE2_MUL( a[4] , b[4] ) , SSE2_MUL( a[5] , b[7] ) ) +
    _mm_add_pd( SSE2_MUL( a[6] , b[2] ) , SSE2_MUL( a[7] , b[5] ) ) +
    SSE2_MUL( a[8] , b[8] ) ;
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
#endif
  return ;
}

// 3x3 mult a = ( b^{\dagger} ).c 
void 
multabdag( __m128d *__restrict a , 
	   const __m128d *__restrict b , 
	   const __m128d *__restrict c )
{
#if NC==3
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
#endif
  return ;
}

// a = b * c^{\dagger}
void 
multab_dag( __m128d *__restrict a , 
	    const __m128d *__restrict b , 
	    const __m128d *__restrict c )
{
#if NC==3
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
#endif
  return ;
}


// a = b^{\dagger} * c^{\dagger}
void 
multab_dagdag( __m128d *__restrict a , 
	       const __m128d *__restrict b , 
	       const __m128d *__restrict c )
{
#if NC==3
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
#endif
  return ;
}

// print a to stdout
void
print_colormatrix( const double complex a[ NCNC ] )
{
  int i , j ;
  printf( "\n" ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      printf( "%f %f " , creal( a[j+i*NC] ) , cimag( a[j+i*NC] ) ) ;
    }
    printf( "\n" ) ;
  }
  return ;
}

#endif
