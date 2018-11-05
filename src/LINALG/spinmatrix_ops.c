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

// computes GLEFT spinmatrix GRIGHT
void
gamma_spinmatrix_lr( struct spinmatrix *S ,
		     const struct gamma GLEFT ,
		     const struct gamma GRIGHT )
{
  struct spinmatrix tmp = *S ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    const size_t col1 = GLEFT.ig[ i ] ;
    for( j = 0 ; j < NS ; j++ ) {
      const size_t col2 = GRIGHT.ig[ j ] ;
      // switch for the phases
      switch( ( GLEFT.g[ i ] + GRIGHT.g[ col2 ] ) & 3 ) {
      case 0 : tmp.D[i][j] = S -> D[col1][col2] ; break ;
      case 1 : tmp.D[i][j] = I * S -> D[col1][col2] ; break ;
      case 2 : tmp.D[i][j] = -S -> D[col1][col2] ; break ;
      case 3 : tmp.D[i][j] = -I * S -> D[col1][col2] ; break ;
      }
    }
  }
  *S = tmp ;
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
  const double complex *dS = (const double complex*)( S.D ) ; 
  double complex *s = (double complex*)spinmatrix ;
  dS += c2 + NC * c1  ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    *s = *dS ; s++ ; dS += NCNC ;
  }
  return ;
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

// computes a = b * c^T
void
spinmatrix_multiply_T( void *a ,
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
	REC = creal( C[ m + NS*j ] ) ; IMC = cimag( C[ m + NS*j ] ) ;
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
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    sum += s[ i * ( NS + 1 ) ] ;
  }
  return sum ;
}

// trace of the product of two spinmatrices
double complex
trace_prod_spinmatrices( const void *a , 
			 const void *b )
{
  const double complex *A = (const double complex*)a ;
  const double complex *B = (const double complex*)b ;
#if NS == 4
  return 
    A[0]  * B[0] + A[1]  * B[4] + A[2]  * B[8]  + A[3]  * B[12] + \
    A[4]  * B[1] + A[5]  * B[5] + A[6]  * B[9]  + A[7]  * B[13] + \
    A[8]  * B[2] + A[9]  * B[6] + A[10] * B[10] + A[11] * B[14] + \
    A[12] * B[3] + A[13] * B[7] + A[14] * B[11] + A[15] * B[15] ;
#else
  size_t i , j ;
  register double complex sum = 0.0 ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum += A[ j + i * NS ] * B[ i + j * NS ] ;
    }
  }
  return sum ;
#endif
}

// trace of the product of two spinmatrices
double complex
trace_prod_spinmatrices_dag( const void *a , 
			     const void *b )
{
  const double complex *A = (const double complex*)a ;
  const double complex *B = (const double complex*)b ;
#if NS == 4
  return 
    A[0]  * B[0] + A[1]  * conj(B[1] + A[2]  * B[8]  + A[3]  * B[12] +	\
    A[4]  * B[1] + A[5]  * B[5] + A[6]  * B[9]  + A[7]  * B[13] + \
    A[8]  * B[2] + A[9]  * B[6] + A[10] * B[10] + A[11] * B[14] + \
    A[12] * B[3] + A[13] * B[7] + A[14] * B[11] + A[15] * B[15] ;
#else
  size_t i , j ;
  register double complex sum = 0.0 ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      sum += A[ j + i * NS ] * B[ i + j * NS ] ;
    }
  }
  return sum ;
#endif
}

// transpose a matrix
void
transpose_spinmatrix( void *a )
{
  double complex *A = (double complex*)a ;
  register double complex tmp ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = i+1 ; j < NS ; j++ ) {
      tmp = A[ j + NS*i ] ;
      A[ j + NS*i ] = A[ i + NS*j ] ;
      A[ i + NS*j ] = tmp ;
    }
  }
}

// zero a spinmatrix
void
zero_spinmatrix( void *a )
{
  double complex *A = (double complex*)a ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *A = 0.0 , A++ ;
  }
  return ;
}

#endif
