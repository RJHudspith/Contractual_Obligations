/**
   @file matrix_ops_SSE.c
   @brief simple matrix multiplies (SSEd version)

   LU decomposition taken from my code GLU - J
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
  *a = _mm_add_pd( *a , *b ) ; 
#elif NC == 2
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; a++ ; b++ ;
  *a = _mm_add_pd( *a , *b ) ; 
#else
  size_t i ;
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
  memcpy( a , b , NCNC*sizeof( double complex ) ) ;
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
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = (double complex)b[i] ;
  }
#endif
  return ;
}

// double matrix to float
void
colormatrix_equiv_d2f( float complex a[ NCNC ] ,
		       const double complex b[ NCNC ] )
{
#if NC == 3
  a[0] = (float complex)b[0] ; a[1] = (float complex)b[1] ; a[2] = (float complex)b[2] ;
  a[3] = (float complex)b[3] ; a[4] = (float complex)b[4] ; a[5] = (float complex)b[5] ;
  a[6] = (float complex)b[6] ; a[7] = (float complex)b[7] ; a[8] = (float complex)b[8] ;
#elif NC == 2
  a[0] = (float complex)b[0] ; a[1] = (float complex)b[1] ;
  a[2] = (float complex)b[2] ; a[3] = (float complex)b[3] ;
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = (float complex)b[i] ;
  }
#endif
  return ;
}

// computes a[i] = S*( a[i] - b[i] )
void
colormatrix_Sa_xmy( double complex a[ NCNC ] ,
		    const double complex b[ NCNC ] ,
		    const double S )
{
  __m128d *pA = (__m128d*)a ;
  const __m128d *pB = (const __m128d*)b ;
  const __m128d s = _mm_set_pd( S , S ) ;
  // this is basically an FMA - should update the instruction
#if NC == 3
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ;
#else
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    *pA =  _mm_mul_pd( s , _mm_sub_pd( *pA , *pB ) ) ; pA++ ; pB++ ;
  }
#endif
}

// computes A[i] = S*(i*b) + A[i]
void
colormatrix_iSaxpy( double complex a[ NCNC ] ,
		   const double complex b[ NCNC ] ,
		   const double S )
{
  __m128d *pA = (__m128d*)a ;
  const __m128d *pB = (const __m128d*)b ;
  const __m128d s = _mm_set_pd( S , S ) ;
#if NC == 3
  #ifdef __FMA__
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ;
  #else
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ;
  #endif
#else
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    #ifdef __FMA__
    *pA = _mm_fmadd_pd( s , SSE2_iMUL( *pB ) , *pA ) ; pA++ ; pB++ ;
    #else
    *pA = _mm_add_pd( _mm_mul_pd( s , SSE2_iMUL( *pB ) ) , *pA ) ; pA++ ; pB++ ;
    #endif
  }
#endif
}

// computes A[i] = S*b + A[i]
void
colormatrix_Saxpy( double complex a[ NCNC ] ,
		   const double complex b[ NCNC ] ,
		   const double S )
{
  __m128d *pA = (__m128d*)a ;
  const __m128d *pB = (const __m128d*)b ;
  const __m128d s = _mm_set_pd( S , S ) ;
#if NC == 3
  #ifdef __FMA__
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_fmadd_pd( s , *pB , *pA ) ;
  #else
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
  *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ;
  #endif
#else
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    #ifdef __FMA__
    *pA = _mm_fmadd_pd( s , *pB , *pA ) ; pA++ ; pB++ ;
    #else
    *pA = _mm_add_pd( _mm_mul_pd( s , *pB ) , *pA ) ; pA++ ; pB++ ;
    #endif
  }
#endif
}

// trace of a color matrix
__m128d
colortrace( const __m128d *a )
{
#if NC == 3
  return _mm_add_pd( _mm_add_pd( a[ 0 ] , a[ 4 ] ) , a[8] ) ; 
#else
  register __m128d sum = _mm_setzero_pd( ) ;
  size_t i ;
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
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
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
  const __m128d *pU = (const __m128d*)U ;
  __m128d *pR = (__m128d*)res ;
  register const __m128d C = _mm_setr_pd( creal( constant ) ,
					  cimag( constant ) ) ;
#if NC == 3
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
  *pR = SSE2_MUL( C , *pU ) ;
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *pR = SSE2_MUL( C , *pU ) ; pU++ ; pR++ ;
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
  *res = SSE2_CONJ( *( U + 8 ) ) ; 
#elif NC == 2
  *res = SSE2_CONJ( *( U + 0 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 2 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 1 ) ) ; res++ ; 
  *res = SSE2_CONJ( *( U + 3 ) ) ; 
#else
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      *res = SSE2_CONJ( U[ i + j * NC ] ) ; res++ ;
    }
  }
#endif
  return ;
}

//  LU determinant computation with intrinsics
double complex
LU_det( const size_t N , const double complex U[ N*N ] )
{
  double complex a[ N*N ] __attribute__((aligned(ALIGNMENT))) ;
  memcpy( a , U , N*N * sizeof( double complex) ) ;
  return LU_det_overwrite( N , a ) ;
}

//  same as above overwrites U
double complex
LU_det_overwrite( const size_t N , double complex U[ N*N ] )
{
  __m128d determinant = _mm_setr_pd( 1.0 , 0.0 ) ;
  register __m128d best , attempt , z1 , z2 ; 
  double complex s ;
  size_t i , j , l , piv , perms = 0 ;
  __m128d *a = (__m128d*)U ;

  for( i = 0 ; i < N-1 ; i++ ) {
    // z1 = a[i*(N+1)] , z2 = Im(z1),Re(z1)
    z1 = a[i*(N+1)] ;
    z2 = _mm_shuffle_pd( z1 , z1 , 1 ) ;
    best = _mm_add_pd( _mm_mul_pd( z1 , z1 ) ,
		       _mm_mul_pd( z2 , z2 ) ) ;
    piv = i ;
    // compare frob norms of other elements
    for( j = i+1 ; j < N ; j++ ) {
      z1 = a[i+j*N] ;
      z2 = _mm_shuffle_pd( z1 , z1 , 1 ) ;
      attempt = _mm_add_pd( _mm_mul_pd( z1 , z1 ) ,
			    _mm_mul_pd( z2 , z2 ) ) ;
      if( _mm_ucomilt_sd( best , attempt ) ) { 
	piv = j ; 
	best = attempt ; 
      }
    }
    if( piv != i ) {
      // physically swap rows
      __m128d *p1 = a + i*N , *p2 = a + piv * N ;
      for( l = 0 ; l < N ; l++ ) {
	z1 = *p1 ; *p1++ = *p2 ; *p2++ = z1 ;
      }
      perms++ ;
    }
    if( _mm_ucomile_sd( best , _mm_setzero_pd() ) ) {
      fprintf( stderr , "[DETERMINANT] LU  Singular Matrix!!!\n" ) ;
      return 0.0 ;
    }
    // perform gaussian elimination
    const __m128d dt = _mm_div_pd( SSE2_CONJ( a[ i*(N+1) ] ) ,
				   best ) ;

    for( j = N-1 ; j > i ; j-- ) { // go up in other column
      __m128d *pB = a + i + j*N ;
      register const __m128d fac1 = SSE2_MUL( *pB , dt ) ; pB++ ;
      // go along the row performing the subtraction, there is no point in
      // subtracting elements where we have determined the best pivot, just the
      // columns to the right of the pivot
      const __m128d *pA = a + i*(N+1) + 1 ;
      for( l = 0 ; l < N - i - 1 ; l++ ) {
	*pB = _mm_sub_pd( *pB , SSE2_MUL( fac1 , *pA ) ) , pB++ , pA++ ;
      }
    }
    determinant = SSE2_MUL( determinant , a[ i*(N+1) ] ) ;
  }
  determinant = SSE2_MUL( determinant , a[ N*N-1 ] ) ;
  _mm_store_pd( (void*)&s , determinant ) ;
  return perms&1 ? -s : s ;
}

// print "a" to stdout
void
print_colormatrix( const double complex a[ NCNC ] )
{
  size_t i , j ;
  fprintf( stdout , "\n" ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      fprintf( stdout , "%f %f " , creal( a[j+i*NC] ) , cimag( a[j+i*NC] ) ) ;
    }
    fprintf( stdout , "\n" ) ;
  }
  return ;
}

// zero a colormatrix
void
zero_colormatrix( const double complex a[ NCNC ] )
{
  __m128d *pA = (__m128d*)a ;
#if NC == 3
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; 
#elif NC == 2
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; pA++ ;
  *pA = _mm_setzero_pd() ; 
#else
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    *pA = _mm_setzero_pd() ; pA++ ;
  }
#endif
}

#endif
