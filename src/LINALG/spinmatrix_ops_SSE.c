/**
   @file spinmatrix_ops_SSE.c
   @brief SSEd spinmatrix operations
 */
#include "common.h"

#ifdef HAVE_EMMINTRIN_H

#include "spinmatrix_ops.h" // include itself for alphabetising

// atomically add two spinmatrices
void
atomic_add_spinmatrices( void *res ,
			 const void *D )
{
  __m128d *r = (__m128d*)res ;
  const __m128d *d = (const __m128d*)D ;
#if NS == 4
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    *r = _mm_add_pd( *r , *d ) ; r++ ; d++ ;
    *r = _mm_add_pd( *r , *d ) ; r++ ; d++ ;
    *r = _mm_add_pd( *r , *d ) ; r++ ; d++ ;
    *r = _mm_add_pd( *r , *d ) ; r++ ; d++ ;
  }
#else
  int i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *r = _mm_add_pd( *r , *d ) ; r++ ; d++ ;
  }
#endif
  return ;
}

// compute p-slash spinmatrix
void
compute_pslash( void *pslash , 
		const struct gamma *GAMMA ,
		const size_t DIMS , 
		const double p[ DIMS ] )
{
  __m128d tmp1[ NSNS ] ;
  __m128d *ps = (__m128d *)pslash ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *ps = _mm_setzero_pd( ) , p++ ;
  }
  // perform the sum
  size_t mu ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    size_t d , j ;
    for( d = 0 ; d < NS ; d++ ) {
      for( j = 0 ; j < NS ; j++ ) {
	tmp1[ j + d*NS ] = ( d == j ) ? p[mu] : 0.0 ; 
      }
    }
    // multiply and add
    gamma_spinmatrix( tmp1 , GAMMA[ mu ] ) ;
    atomic_add_spinmatrices( ps , tmp1 ) ;
  }
  return ;
}

// left multiply a spinmatrix by a gamma G
void
gamma_spinmatrix( void *spinmatrix ,
		  const struct gamma G ) 
{
  __m128d r[ NSNS ] ;
  const __m128d *d = (const __m128d*)spinmatrix ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    const uint8_t col = G.ig[i] ;
    switch( G.g[i] ) {
    case 0 : 
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = d[ j + col * NS ] ;
      }
      break ;
    case 1 : 
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE2_iMUL( d[ j + col * NS ] ) ;
      }
      break ;
    case 2 : 
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE_FLIP( d[ j + col * NS ] ) ;
      }
      break ;
    case 3 : 
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE2_miMUL( d[ j + col * NS ] ) ;
      }
      break ;
    }
  }
  memcpy( (void*)d , r , NSNS * sizeof( __m128d ) ) ;
  return ;
}

// computes spintrace :: Tr[ G spinmatrix ]
double complex
gammaspinmatrix_trace( const struct gamma G ,
		       const void *spinmatrix )
{
  const __m128d *D = (const __m128d*)spinmatrix ;
  register __m128d sum = _mm_setzero_pd( ) ;
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    switch( G.g[i] ) {
    case 0 : 
      sum = _mm_add_pd( sum , D[ i + G.ig[i] * NS ] ) ; 
      break ;
    case 1 : 
      sum = _mm_add_pd( sum , SSE2_iMUL( D[ i + G.ig[i] * NS ] ) ) ; 
      break ;
    case 2 : 
      sum = _mm_add_pd( sum , SSE_FLIP( D[ i + G.ig[i] * NS ] ) ) ; 
      break ;
    case 3 : 
      sum = _mm_add_pd( sum , SSE2_miMUL( D[ i + G.ig[i] * NS ] ) ) ; 
      break ;
    }
  }
  double complex csum ;
  _mm_store_pd( (void*)&csum , sum ) ;
  return csum ;
}

// right multiply a spinmatrix by a gamma G
void
spinmatrix_gamma( void *spinmatrix ,
		  const struct gamma G ) 
{
  __m128d r[ NSNS ] ; // temporary space
  const __m128d *d = (const __m128d*)spinmatrix ;
  size_t i , j ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const uint8_t col = G.ig[j] ;
    switch( G.g[ col ] ) {
    case 0 : 
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = d[ col + i * NS ] ;
      } 
      break ;
    case 1 : 
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE2_iMUL( d[ col + i * NS ] ) ;
      }
      break ;
    case 2 : 
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE_FLIP( d[ col + i * NS ] ) ;
      }
      break ;
    case 3 :
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE2_miMUL( d[ col + i * NS ] ) ;
      }
      break ;
    }
  }
  memcpy( (void*)d , r , NSNS * sizeof( __m128d ) ) ;
  return ;
}

// computes a = b * c
void
spinmatrix_multiply( void *a ,
		     const void *b ,
		     const void *c )
{
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
  size_t i , j ;
  register __m128d sum ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      #if NS == 4
      sum = _mm_add_pd( sum , SSE2_MUL( B[ 0 + i * NS ] , C[ j ] ) );
      sum = _mm_add_pd( sum , SSE2_MUL( B[ 1 + i * NS ] , C[ j + 1 * NS ] ) );
      sum = _mm_add_pd( sum , SSE2_MUL( B[ 2 + i * NS ] , C[ j + 2 * NS ] ) );
      sum = _mm_add_pd( sum , SSE2_MUL( B[ 3 + i * NS ] , C[ j + 3 * NS ] ) );
      #else
      size_t m ;
      for( m = 0 ; m < NS ; m++  ) {
	sum = _mm_add_pd( sum , SSE2_MUL( B[ m + i * NS ] , C[ j + m * NS ] ) );
      }
      #endif
      *A = sum , A++ ;
    }
  }
  return ;
}

// atomically multiply a spinmatrix by a constant factor
void
spinmatrix_mulconst( void *spinmatrix , 
		     const double factor )
{
  __m128d *s = (__m128d*)spinmatrix ;
  register const __m128d f = _mm_set_pd( factor , factor ) ;
  size_t i ;
#if NS == 4
  for( i = 0 ; i < NS ; i++ ) {
    *s = _mm_mul_pd( *s , f ) ; s++ ;
    *s = _mm_mul_pd( *s , f ) ; s++ ;
    *s = _mm_mul_pd( *s , f ) ; s++ ;
    *s = _mm_mul_pd( *s , f ) ; s++ ;
  }
#else
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = _mm_mul_pd( *s , f ) ; s++ ;
  }
#endif
  return ;
}

// trace of a spinmatrix
double complex
spinmatrix_trace( const void *spinmatrix )
{
  const __m128d *s = (const __m128d*)spinmatrix ;
  register __m128d sum = _mm_setzero_pd( ) ;
#if NS == 4
  sum = _mm_add_pd( sum , *s ) , s+= NS + 1 ;
  sum = _mm_add_pd( sum , *s ) , s+= NS + 1 ;
  sum = _mm_add_pd( sum , *s ) , s+= NS + 1 ;
  sum = _mm_add_pd( sum , *s ) , s+= NS + 1 ;
#else
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    sum += s[ i * ( NS + 1 ) ] ;
  }
#endif
  double complex csum ;
  _mm_store_pd( (void*)&csum , sum ) ;
  return csum ;
}

#endif
