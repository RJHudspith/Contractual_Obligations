/**
   @file spinmatrix_ops.c
   @brief operations acting on #NSNS flattened spin-matrices
 */
#include "common.h"

#ifndef HAVE_EMMINTRIN_H

#include "spinmatrix_ops.h" // include itself

// atomically add two spinmatrices
void
atomic_add_spinmatrices( void *res ,
			 const void *D )
{
  double complex *r = (double complex*)res ;
  const double complex *d = (const double complex*)D ;
  int i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    r[ i ] += d[ i ] ;
  }
  return ;
}

// compute p-slash spinmatrix
void
compute_pslash( void *pslash , 
		const struct gamma *GAMMA ,
		const double p[ ND ] )
{
  double complex *ps = (double complex*)pslash ;
  size_t mu , s ;
  // set pslash to 0.0
  for( mu = 0 ; mu < NSNS ; mu++ ) {
    ps[ mu ] = 0.0 ;
  }
  // perform the sum
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( s = 0 ; s < NS ; s++ ) {
      switch( GAMMA[ mu ].g[ s ] ) {
      case 0 : ps[ GAMMA[ mu ].ig[ s ] + s * NS ] +=      p[mu] ; break ;
      case 1 : ps[ GAMMA[ mu ].ig[ s ] + s * NS ] +=  I * p[mu] ; break ;
      case 2 : ps[ GAMMA[ mu ].ig[ s ] + s * NS ] -=      p[mu] ; break ;
      case 3 : ps[ GAMMA[ mu ].ig[ s ] + s * NS ] += -I * p[mu] ; break ;
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
  //double complex *r = (double complex*)res ;
  double complex r[ NSNS ] ;
  double complex *d = (double complex*)spinmatrix ;
  int i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    const int col = (int)G.ig[i] ;
    register double complex factor = 1.0 ;
    switch( G.g[i] ) {
    case 0 : factor = +1 ; break ;
    case 1 : factor = +I ; break ;
    case 2 : factor = -1 ; break ;
    case 3 : factor = -I ; break ;
    }
    for( j = 0 ; j < NS ; j++ ) {
      r[ j + i * NS ] = factor * d[ j + col * NS ] ;
    }
  }
  memcpy( d , r , NSNS * sizeof( double complex ) ) ;
  return ;
}

// computes spintrace :: Tr[ G spinmatrix ]
double complex
gammaspinmatrix_trace( const struct gamma G ,
		       const void *spinmatrix )
{
  const double complex *D = (const double complex*)spinmatrix ;
  register double complex sum = 0.0 ;
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    switch( G.g[i] ) {
    case 0 : sum += D[ i + G.ig[i] * NS ] ; break ;
    case 1 : sum += +I * D[ i + G.ig[i] * NS ] ; break ;
    case 2 : sum += -D[ i + G.ig[i] * NS ] ; break ;
    case 3 : sum += -I * D[ i + G.ig[i] * NS ] ; break ;
    }
  }
  return sum ;
}

// set the spinmatrix to the identity
void
identity_spinmatrix( void *spinmatrix )
{
  double complex *s = (double complex*)spinmatrix ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      *s = ( i != j ) ? 0.0 : 1.0 ; s++ ;
    }
  }
  return ;
}

// right multiply a spinmatrix by a gamma G
void
spinmatrix_gamma( void *spinmatrix ,
		  const struct gamma G ) 
{
  double complex r[ NSNS ] ; // temporary space
  double complex *d = (double complex*)spinmatrix ;
  size_t i , j ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const int col = G.ig[j] ;
    double complex fac = 1.0 ;
    switch( G.g[ col ] ) {
    case 0 : fac =  1 ; break ;
    case 1 : fac =  I ; break ;
    case 2 : fac = -1 ; break ;
    case 3 : fac = -I ; break ;
    }
    // and multiply it in
    for( i = 0 ; i < NS ; i++ ) {
      r[ j + i * NS ] = fac * d[ col + i * NS ] ;
    } 
  }
  memcpy( d , r , NSNS * sizeof( double complex ) ) ;
  return ;
}

// computes a = b * c
void
spinmatrix_multiply( void *a ,
		     const void *b ,
		     const void *c )
{
  double complex *A = (double complex*)a ;
  const double complex *B = (double complex*)b ;
  const double complex *C = (double complex*)c ;
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NS ; m++  ) {
	REB = creal( B[ m + NS*i ] ) ; IMB = cimag( B[ m + NS*i ] ) ;
	REC = creal( C[ j + m*NS ] ) ; IMC = cimag( C[ j + m*NS ] ) ;
	sum += REB * REC - IMB * IMC + I * ( REB * IMC + IMB * REC ) ; 
      }
      A[ j + NS*i ] = sum ;
    }
  }
  return ;
}

// atomically multiply a spinmatrix by a constant factor
void
spinmatrix_mulconst( void *spinmatrix , 
		     const double factor )
{
  double complex *s = (double complex*)spinmatrix ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    s[ i ] *= factor ;
  }
  return ;
}

// trace of a spinmatrix
double complex
spinmatrix_trace( const void *spinmatrix )
{
  const double complex *s = (const double complex *)spinmatrix ;
  register double complex sum = 0.0 ;
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    sum += s[ i * ( NS + 1 ) ] ;
  }
  return sum ;
}

// zero a spinmatrix
void
zero_spinmatrix( void *spinmatrix )
{
  double complex *s = (double complex*)spinmatrix ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *s = 0.0 , s++ ;
  }
  return ;
}

#endif
