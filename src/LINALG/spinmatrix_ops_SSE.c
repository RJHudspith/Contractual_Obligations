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
		const double p[ ND ] )
{
  __m128d *ps = (__m128d *)pslash ;
  size_t mu , s ;
  // set pslash to 0.0
  for( mu = 0 ; mu < NSNS ; mu++ ) {
    ps[ mu ] = _mm_setzero_pd( ) ;
  }
  // perform the sum
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( s = 0 ; s < NS ; s++ ) {
      switch( GAMMA[ mu ].g[ s ] ) {
      case 0 : 
	ps[ GAMMA[ mu ].ig[ s ] + s * NS ] = 
	  _mm_add_pd( ps[ GAMMA[ mu ].ig[ s ] + s * NS ] , 
		      _mm_setr_pd( p[mu] , 0 ) ) ; 
	break ;
      case 1 : 
	ps[ GAMMA[ mu ].ig[ s ] + s * NS ] = 
	  _mm_add_pd( ps[ GAMMA[ mu ].ig[ s ] + s * NS ] , 
		      _mm_setr_pd( 0 , p[mu] ) ) ; 
	break ;
      case 2 : 
	ps[ GAMMA[ mu ].ig[ s ] + s * NS ] = 
	  _mm_add_pd( ps[ GAMMA[ mu ].ig[ s ] + s * NS ] , 
		      _mm_setr_pd( -p[mu] , 0 ) ) ; 
	break ;
      case 3 : 
	ps[ GAMMA[ mu ].ig[ s ] + s * NS ] = 
	  _mm_add_pd( ps[ GAMMA[ mu ].ig[ s ] + s * NS ] , 
		      _mm_setr_pd( 0 , -p[mu] ) ) ; 
	break ;
      }
    }
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
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    const uint8_t col = G.ig[i] ;
    switch( G.g[i] ) {
    case 0 :
      #if NS == 4
      r[ 0 + i * NS ] = d[ 0 + col * NS ] ;
      r[ 1 + i * NS ] = d[ 1 + col * NS ] ;
      r[ 2 + i * NS ] = d[ 2 + col * NS ] ;
      r[ 3 + i * NS ] = d[ 3 + col * NS ] ;
      #else
      size_t j ;
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = d[ j + col * NS ] ;
      }
      #endif
      break ;
    case 1 :
      #if NS == 4
      r[ 0 + i * NS ] = SSE2_iMUL( d[ 0 + col * NS ] ) ;
      r[ 1 + i * NS ] = SSE2_iMUL( d[ 1 + col * NS ] ) ;
      r[ 2 + i * NS ] = SSE2_iMUL( d[ 2 + col * NS ] ) ;
      r[ 3 + i * NS ] = SSE2_iMUL( d[ 3 + col * NS ] ) ;
      #else
      size_t j ;
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE2_iMUL( d[ j + col * NS ] ) ;
      }
      #endif
      break ;
    case 2 : 
      #if NS == 4
      r[ 0 + i * NS ] = SSE_FLIP( d[ 0 + col * NS ] ) ;
      r[ 1 + i * NS ] = SSE_FLIP( d[ 1 + col * NS ] ) ;
      r[ 2 + i * NS ] = SSE_FLIP( d[ 2 + col * NS ] ) ;
      r[ 3 + i * NS ] = SSE_FLIP( d[ 3 + col * NS ] ) ;
      #else
      size_t j ;
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE_FLIP( d[ j + col * NS ] ) ;
      }
      #endif
      break ;
    case 3 : 
      #if NS == 4
      r[ 0 + i * NS ] = SSE2_miMUL( d[ 0 + col * NS ] ) ;
      r[ 1 + i * NS ] = SSE2_miMUL( d[ 1 + col * NS ] ) ;
      r[ 2 + i * NS ] = SSE2_miMUL( d[ 2 + col * NS ] ) ;
      r[ 3 + i * NS ] = SSE2_miMUL( d[ 3 + col * NS ] ) ;
      #else
      size_t j ;
      for( j = 0 ; j < NS ; j++ ) {
	r[ j + i * NS ] = SSE2_miMUL( d[ j + col * NS ] ) ;
      }
      #endif
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

// gets a spinmatrix from our propagator
void
get_spinmatrix( void *spinmatrix , 
		const struct spinor S ,
		const size_t c1 ,
		const size_t c2 )
{
  // pokes specific colors c1,c2 out from the spinor
  // each spinor is 4x4x3x3 chooses a component of the 3x3
  // color matrix and hence gives a 4x4 spinmatrix
  const __m128d *dS = (const __m128d*)( S.D ) ; 
  __m128d *s = (__m128d*)spinmatrix ;
  dS += c2 + NC * c1  ; // jump to the first colour component
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    *s = *dS ; s++ ; 
    dS += NCNC ; // jump pointer forward a full colour matrix-worth of data
  }
  return ;
}

// set spinmatrix to the identity
void 
identity_spinmatrix( void *spinmatrix )
{
  __m128d *s = (__m128d*)spinmatrix ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      *s = ( i != j ) ? _mm_setzero_pd() : \
	_mm_setr_pd( 1.0 , 0.0 ) ; s++ ;
    }
  }
  return ;
}

// right multiply a spinmatrix by a gamma G
void
spinmatrix_gamma( void *spinmatrix ,
		  const struct gamma G ) 
{
  __m128d r[ NSNS ] ; // temporary space
  const __m128d *d = (const __m128d*)spinmatrix ;
  size_t j ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const uint8_t col = G.ig[j] ;
    switch( G.g[ col ] ) {
    case 0 : 
      #if NS == 4
      r[ j + 0  ] = d[ col + 0  ] ;
      r[ j + 4  ] = d[ col + 4  ] ;
      r[ j + 8  ] = d[ col + 8  ] ;
      r[ j + 12 ] = d[ col + 12 ] ;
      #else
      size_t i ;
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = d[ col + i * NS ] ;
      } 
      #endif
      break ;
    case 1 : 
      #if NS == 4
      r[ j + 0  ] = SSE2_iMUL( d[ col + 0  ] ) ;
      r[ j + 4  ] = SSE2_iMUL( d[ col + 4  ] ) ;
      r[ j + 8  ] = SSE2_iMUL( d[ col + 8  ] ) ;
      r[ j + 12 ] = SSE2_iMUL( d[ col + 12 ] ) ;
      #else
      size_t i ;
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE2_iMUL( d[ col + i * NS ] ) ;
      }
      #endif
      break ;
    case 2 :
      #if NS == 4
      r[ j + 0  ] = SSE_FLIP( d[ col + 0  ] ) ;
      r[ j + 4  ] = SSE_FLIP( d[ col + 4  ] ) ;
      r[ j + 8  ] = SSE_FLIP( d[ col + 8  ] ) ;
      r[ j + 12 ] = SSE_FLIP( d[ col + 12 ] ) ;
      #else
      size_t i ;
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE_FLIP( d[ col + i * NS ] ) ;
      }
      #endif
      break ;
    case 3 :
      #if NS == 4
      r[ j + 0  ] = SSE2_miMUL( d[ col + 0  ] ) ;
      r[ j + 4  ] = SSE2_miMUL( d[ col + 4  ] ) ;
      r[ j + 8  ] = SSE2_miMUL( d[ col + 8  ] ) ;
      r[ j + 12 ] = SSE2_miMUL( d[ col + 12 ] ) ;
      #else
      size_t i ;
      for( i = 0 ; i < NS ; i++ ) {
	r[ j + i * NS ] = SSE2_miMUL( d[ col + i * NS ] ) ;
      }
      #endif
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
#if NS == 4
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x4) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xc) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x1) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0x9) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xd) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x6) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xe) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x3) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x7) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0xb) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xf) ) ) ) ;	    
  A++ ;
  // next row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x4) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xc) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x1) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0x9) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xd) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x6) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xe) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x3) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x7) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0xb) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xf) ) ) ) ;
  A++ ;
  // third row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x4) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xc) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x1) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0x9) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xd) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x6) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xe) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x3) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x7) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0xb) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xf) ) ) ) ;
  A++ ;
  // fourth row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x4) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xc) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x1) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0x9) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xd) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x6) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xe) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x3) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x7) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0xb) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xf) ) ) ) ;
#else
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
  register __m128d sum ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      size_t m ;
      for( m = 0 ; m < NS ; m++  ) {
	sum = _mm_add_pd( sum , SSE2_MUL( B[ m + i * NS ] , C[ j + m * NS ] ) );
      }
      *A = sum , A++ ;
    }
  }
#endif
  return ;
}

// computes a = b * c^T
extern inline void
spinmatrix_multiply_T( void *a ,
		       const void *b ,
		       const void *c )
{
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
#if NS == 4
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x1) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0x3) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x4) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0x6) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0x7) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0x9) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xb) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x0) , *(C + 0xc) ) ,
			       SSE2_MUL( *(B + 0x1) , *(C + 0xd) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x2) , *(C + 0xe) ) ,
			       SSE2_MUL( *(B + 0x3) , *(C + 0xf) ) ) ) ;	    
  A++ ;
  // next row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x1) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0x3) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x4) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0x6) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0x7) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0x9) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xb) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x4) , *(C + 0xc) ) ,
			       SSE2_MUL( *(B + 0x5) , *(C + 0xd) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0x6) , *(C + 0xe) ) ,
			       SSE2_MUL( *(B + 0x7) , *(C + 0xf) ) ) ) ;
  A++ ;
  // third row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x1) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0x3) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x4) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0x6) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0x7) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0x9) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xb) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0x8) , *(C + 0xc) ) ,
			       SSE2_MUL( *(B + 0x9) , *(C + 0xd) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xa) , *(C + 0xe) ) ,
			       SSE2_MUL( *(B + 0xb) , *(C + 0xf) ) ) ) ;
  A++ ;
  // fourth row
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x0) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x1) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0x2) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0x3) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x4) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0x6) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0x7) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0x8) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0x9) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0xa) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xb) ) ) ) ;
  A++ ;
  *A = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(B + 0xc) , *(C + 0xc) ) ,
			       SSE2_MUL( *(B + 0xd) , *(C + 0xd) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(B + 0xe) , *(C + 0xe) ) ,
			       SSE2_MUL( *(B + 0xf) , *(C + 0xf) ) ) ) ;
#else
  __m128d *A = (__m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;
  register __m128d sum ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum = _mm_setzero_pd( ) ;
      size_t m ;
      for( m = 0 ; m < NS ; m++  ) {
	sum = _mm_add_pd( sum , SSE2_MUL( B[ m + i * NS ] , C[ m + j * NS ] ) );
      }
      *A = sum , A++ ;
    }
  }
#endif
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
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    sum += s[ i * ( NS + 1 ) ] ;
  }
#endif
  double complex csum ;
  _mm_store_pd( (void*)&csum , sum ) ;
  return csum ;
}

// trace of the product of two spinmatrices
double complex
trace_prod_spinmatrices( const void *a , 
			 const void *b )
{
  const __m128d *A = (const __m128d*)a ;
  const __m128d *B = (const __m128d*)b ;
  double complex csum ;
#if NS == 4
  register __m128d s1 , s2 , s3 , s4 ;
  s1 = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(A + 0x0) , *(B + 0x0) ) , 
			       SSE2_MUL( *(A + 0x1) , *(B + 0x4) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(A + 0x2) , *(B + 0x8) ) , 
			       SSE2_MUL( *(A + 0x3) , *(B + 0xc) ) ) ) ;
  s2 = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(A + 0x4) , *(B + 0x1) ) , 
			       SSE2_MUL( *(A + 0x5) , *(B + 0x5) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(A + 0x6) , *(B + 0x9) ) , 
			       SSE2_MUL( *(A + 0x7) , *(B + 0xd) ) ) ) ;
  s3 = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(A + 0x8) , *(B + 0x2) ) , 
			       SSE2_MUL( *(A + 0x9) , *(B + 0x6) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(A + 0xa) , *(B + 0xa) ) , 
			       SSE2_MUL( *(A + 0xb) , *(B + 0xe) ) ) ) ;
  s4 = _mm_add_pd( _mm_add_pd( SSE2_MUL( *(A + 0xc) , *(B + 0x3) ) , 
			       SSE2_MUL( *(A + 0xd) , *(B + 0x7) ) ) ,
		   _mm_add_pd( SSE2_MUL( *(A + 0xe) , *(B + 0xb) ) , 
			       SSE2_MUL( *(A + 0xf) , *(B + 0xf) ) ) ) ;
  // put them all in s1
  s1 = _mm_add_pd( _mm_add_pd( s1 , s2 ) , _mm_add_pd( s3 , s4 ) ) ;
  _mm_store_pd( (void*)&csum , s1 ) ;
#else
  size_t i , j ;
  register double complex sum = _mm_setzero_pd( ) ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum = _mm_add_pd( sum , SSE2_MUL( A[ j + i * NS ] , 
					B[ i + j * NS ] ) ) ;
    }
  }
  _mm_store_pd( (void*)&csum , sum) ;
#endif
  return csum ;
}

// transpose a spinmatrix
void
transpose_spinmatrix( void *a )
{
  __m128d *A = (__m128d*)a ;
  register __m128d tmp ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = i+1 ; j < NS ; j++ ) {
      tmp = A[ j + NS*i ] ;
      A[ j + NS*i ] = A[ i + NS*j ] ;
      A[ i + NS*j ] = tmp ;
    }
  }
}

// set spinmatrix to zero
void
zero_spinmatrix( void *a )
{
  __m128d *A = (__m128d*)a ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *A = _mm_setzero_pd( ) , A++ ;
  }
  return ;
}

#endif
