/**
   @file matrix_ops.c
   @brief simple matrix multiplies
 */

#include "common.h"
#include "matrix_ops.h"

#ifndef HAVE_EMMINTRIN_H

// add two color matrices
void
add_mat( double complex *__restrict a ,
	 const double complex *__restrict b )
{
#if NC == 3
  a[0] += b[0] ; a[1] += b[1] ; a[2] += b[2] ; 
  a[3] += b[3] ; a[4] += b[4] ; a[5] += b[5] ; 
  a[6] += b[6] ; a[7] += b[7] ; a[8] += b[8] ; 
#elif NC == 2
  a[0] += b[0] ; a[1] += b[1] ;
  a[2] += b[2] ; a[3] += b[3] ;
#else
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] += b[ i ] ;
  }
#endif
}

// equate two color matrices
void
colormatrix_equiv( double complex *__restrict a ,
		   const double complex *__restrict b )
{
#if NC == 3
  a[0] = b[0] ; a[1] = b[1] ; a[2] = b[2] ; 
  a[3] = b[3] ; a[4] = b[4] ; a[5] = b[5] ; 
  a[6] = b[6] ; a[7] = b[7] ; a[8] = b[8] ;
#elif NC == 2
  a[0] = b[0] ; a[1] = b[1] ;
  a[2] = b[2] ; a[3] = b[3] ;
#else
  memcpy( a , b , NCNC * sizeof( double complex ) ) ;
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

// compute the trace
double complex
colortrace( const double complex a[ NCNC ] )
{
#if NC == 3
  return a[ 0 ] + a[ 4 ] + a[ 8 ] ;
#elif NC == 2
  return a[ 0 ] + a[ 3 ] ;
#else
  register double re = 0. , im = 0. ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    re += creal( a[ i*( NC + 1 ) ] ) ;
    im += cimag( a[ i*( NC + 1 ) ] ) ;
  }
  return re + I * im ;
#endif
}

// is just Tr( a * b )
double complex
colortrace_prod( const double complex *__restrict a , 
		 const double complex *__restrict b )
{
#if NC == 3
  return a[0] * b[0] + a[1] * b[3] + a[2] * b[6] +	\
    a[3] * b[1] + a[4] * b[4] + a[5] * b[7] +		\
    a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;
#elif NC == 2
  return a[0] * b[0] + a[1] * b[2] + a[2] * b[1] + a[3] * b[3] ;
#else
  register GLU_real sumr = 0.0 , sumi = 0.0 ;
  int i, j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sumr += creal( a[ j + NC * i ] ) * creal( b[ i + NC * j ] ) -
	cimag( a[ j + NC * i ] ) * cimag( b[ i + NC * j ] ) ;
      sumi += creal( a[ j + NC * i ] ) * cimag( b[ i + NC * j ] ) +
	cimag( a[ j + NC * i ] ) * creal( b[ i + NC * j ] ) ;
    }
  }
  return sumr + I * sumi ;
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
dagger_gauge( double complex *__restrict res ,
	      const double complex *__restrict U )
{
#if NC == 3
  res[0] = conj( U[0] ) ; res[1] = conj( U[3] ) ; res[2] = conj( U[6] ) ;
  res[3] = conj( U[1] ) ; res[4] = conj( U[4] ) ; res[5] = conj( U[7] ) ;
  res[6] = conj( U[2] ) ; res[7] = conj( U[5] ) ; res[8] = conj( U[8] ) ;
#elif NC == 2
  res[0] = conj( U[0] ) ; res[1] = conj( U[2] ) ; 
  res[2] = conj( U[1] ) ; res[3] = conj( U[3] ) ; 
#else
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      res[ j + i * NC ] = conj( U[ i + j * NC ] ) ;
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
multab( double complex a[ NCNC ] , 
	const double complex b[ NCNC ] , 
	const double complex c[ NCNC ] )
{
#if NC == 3
  a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ;	\
  a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ;	\
  a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ;	\
  a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ;	\
  a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ;	\
  a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ;	\
  a[6] = b[6] * c[0] + b[7] * c[3] + b[8] * c[6] ;	\
  a[7] = b[6] * c[1] + b[7] * c[4] + b[8] * c[7] ;	\
  a[8] = b[6] * c[2] + b[7] * c[5] + b[8] * c[8] ;	
#elif NC == 2
  a[0] = b[0] * c[0] + b[1] * c[2] ;		\
  a[1] = b[0] * c[1] + b[1] * c[3] ;		\
  a[2] = b[2] * c[0] + b[3] * c[2] ;		\
  a[3] = b[2] * c[1] + b[3] * c[3] ;		
#else
  // slow and stupid version
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++  ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ j + m*NC ] ) ; IMC = cimag( c[ j + m*NC ] ) ;
	sum += REB * REC - IMB * IMC + I * ( REB * IMC + IMB * REC ) ; 
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}

// 3x3 mult a = ( b^{\dagger} ).c 
void 
multabdag( double complex a[ NCNC ] , 
	   const double complex b[ NCNC ] , 
	   const double complex c[ NCNC ] )
{
#if NC == 3
  a[0] = conj( b[0] ) * c[0] + conj( b[3] ) * c[3] + conj( b[6] ) * c[6] ; \
  a[1] = conj( b[0] ) * c[1] + conj( b[3] ) * c[4] + conj( b[6] ) * c[7] ; \
  a[2] = conj( b[0] ) * c[2] + conj( b[3] ) * c[5] + conj( b[6] ) * c[8] ; \
  a[3] = conj( b[1] ) * c[0] + conj( b[4] ) * c[3] + conj( b[7] ) * c[6] ; \
  a[4] = conj( b[1] ) * c[1] + conj( b[4] ) * c[4] + conj( b[7] ) * c[7] ; \
  a[5] = conj( b[1] ) * c[2] + conj( b[4] ) * c[5] + conj( b[7] ) * c[8] ; \
  a[6] = conj( b[2] ) * c[0] + conj( b[5] ) * c[3] + conj( b[8] ) * c[6] ; \
  a[7] = conj( b[2] ) * c[1] + conj( b[5] ) * c[4] + conj( b[8] ) * c[7] ; \
  a[8] = conj( b[2] ) * c[2] + conj( b[5] ) * c[5] + conj( b[8] ) * c[8] ;
#elif NC == 2
  a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;	\
  a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;	\
  a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2] ;	\
  a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3] ;
#else
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ i + NC*m ] ) ; IMB = cimag( b[ i + NC*m ] ) ;
	REC = creal( c[ j + NC*m ] ) ; IMC = cimag( c[ j + NC*m ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REB * IMC - IMB * REC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}

// a = b * c^{\dagger}
void 
multab_dag( double complex a[ NCNC ] , 
	    const double complex b[ NCNC ] , 
	    const double complex c[ NCNC ] )
{
#if NC == 3
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; \
  a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; \
  a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; \
  a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; \
  a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; \
  a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ; \
  a[6] = b[6] * conj( c[0] ) + b[7] * conj( c[1] ) + b[8] * conj( c[2] ) ; \
  a[7] = b[6] * conj( c[3] ) + b[7] * conj( c[4] ) + b[8] * conj( c[5] ) ; \
  a[8] = b[6] * conj( c[6] ) + b[7] * conj( c[7] ) + b[8] * conj( c[8] ) ; 
#elif NC == 2
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) ;	\
  a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) ;	\
  a[2] = b[2] * conj( c[0] ) + b[3] * conj( c[1] ) ;	\
  a[3] = b[2] * conj( c[2] ) + b[3] * conj( c[3] ) ;
#else // instead of inlining we have a function call
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0. ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ m + NC*j ] ) ; IMC = cimag( c[ m + NC*j ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REC * IMB - REB * IMC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}


// a = b^{\dagger} * c^{\dagger}
void 
multab_dagdag( double complex a[ NCNC ] , 
	       const double complex b[ NCNC ] , 
	       const double complex c[ NCNC ] )
{
#if NC == 3
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] ) + conj( b[6] ) * conj( c[2] ) ; \
  a[1] = conj( b[0] ) * conj( c[3] ) + conj( b[3] ) * conj( c[4] ) + conj( b[6] ) * conj( c[5] ) ; \
  a[2] = conj( b[0] ) * conj( c[6] ) + conj( b[3] ) * conj( c[7] ) + conj( b[6] ) * conj( c[8] ) ; \
  a[3] = conj( b[1] ) * conj( c[0] ) + conj( b[4] ) * conj( c[1] ) + conj( b[7] ) * conj( c[2] ) ; \
  a[4] = conj( b[1] ) * conj( c[3] ) + conj( b[4] ) * conj( c[4] ) + conj( b[7] ) * conj( c[5] ) ; \
  a[5] = conj( b[1] ) * conj( c[6] ) + conj( b[4] ) * conj( c[7] ) + conj( b[7] ) * conj( c[8] ) ; \
  a[6] = conj( b[2] ) * conj( c[0] ) + conj( b[5] ) * conj( c[1] ) + conj( b[8] ) * conj( c[2] ) ; \
  a[7] = conj( b[2] ) * conj( c[3] ) + conj( b[5] ) * conj( c[4] ) + conj( b[8] ) * conj( c[5] ) ; \
  a[8] = conj( b[2] ) * conj( c[6] ) + conj( b[5] ) * conj( c[7] ) + conj( b[8] ) * conj( c[8] ) ; 
#elif NC == 2
  a[0] = conj( b[0] ) * conj( c[0] ) + conj( b[2] ) * conj( c[1] )  ;	\
  a[1] = conj( b[0] ) * conj( c[2] ) + conj( b[2] ) * conj( c[3] )  ;	\
  a[2] = conj( b[1] ) * conj( c[0] ) + conj( b[3] ) * conj( c[1] )  ;	\
  a[3] = conj( b[1] ) * conj( c[2] ) + conj( b[3] ) * conj( c[3] )  ; 
#else
  int i , j , m ;
  register double sumr = 0.0 , sumi = 0.0 ;
  register double REB , IMB , REC , IMC ;
  const double complex *pC , *pB ;
  for( i = 0 ; i < NC ; i++ ) {
    pC = c ;
    for( j = 0 ; j < NC ; j++ ) {
      pB = b ;
      sumr = sumi = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( pB[i] ) ; IMB = cimag( pB[i] ) ;
	REC = creal( pC[m] ) ; IMC = cimag( pC[m] ) ;
	sumr += REB * REC - IMB * IMC ;
	sumi += REB * IMC + IMB * REC ;
	pB += NC ;
      }
      a[ j + NC*i ] = sumr + I * sumi ;
      pC += NC ;
    }
  }
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