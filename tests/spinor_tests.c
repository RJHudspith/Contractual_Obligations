/**
   @file spinor_tests.c
   @brief spinor_ops testing
 */
#include "common.h"

#include "gammas.h"
#include "minunit.h"     // minimal unit testing framework
#include "spinor_ops.h"  // spinor operations

static struct spinor *S = NULL ; // LCU spinor
static struct spinor A ;         // single spinor storage
static double complex *C ;
static double complex *D ;

#define FTOL ( NC * 1.E-14 ) 

static char
*colortrace_spinor_test( void )
{
  // take the color trace
  colortrace_spinor( C , A.D ) ;
  int j ;
  for( j = 0 ; j < NSNS ; j++ ) {
    // result is an identity
    const size_t res = NC * ( NCNC - 1 ) / 2 + j * NC * ( NCNC ) ;
    mu_assert( "[UNIT] error : spinor ops colortrace_spinor broken " , 
	       !( fabs( creal( C[j] ) - res ) > FTOL 
		  ||
		  fabs( cimag( C[j] ) - res ) > FTOL 
		  ) ) ;
  }
  return NULL ;
}

static char 
*equate_minus_test( void )
{
  struct spinor B ;
  equate_spinor_minus( &B , (const void*)A.D ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor ops equate minus broken " , 
	       !(*a != -*b) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

static char
*flipsign_test( void )
{
  struct spinor B = A ;
  flipsign_spinor( &B ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor ops equate minus broken " , 
	       !(*a != -*b) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

// multiplies spinor by identity matrix and compares
static char
*gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  gauge_spinor( &B , Id , A ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : gauge_spinor broken " , 
	       !( *a != *b ) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

// multiply by (-I)^{\dagger}
static char
*gaugedag_test( void )
{
  double complex mI[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      mI[ j + i * NC ] = ( j != i ) ? 0.0 : -I ;
    }
  }
  struct spinor B ;
  gaugedag_spinor( &B , mI , A ) ;
  const double complex *a = (const double complex*)A.D ;
  const double complex *b = (const double complex*)B.D ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : gaugedag_spinor broken " , 
	       !( fabs( creal( b[i] ) + cimag( a[i] ) ) > FTOL ||
		 fabs( cimag( b[i] ) - creal( a[i] ) ) > FTOL ) ) ;
  }
  return NULL ;
}

// Id * ( A )^{\dagger} where the dagger is JUST OVER COLOR INDICES!
static char
*gauge_spinordag_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  gauge_spinordag( &B , Id , A ) ;
  // B should be the conjugate of the full 12x12 spinor matrix
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      const double complex *a = (const double complex*)A.D[d1][d2].C ;
      const double complex *b = (const double complex*)B.D[d1][d2].C ;
      int c1 , c2 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {			     
	  mu_assert( "[UNIT] error : gauge_spinordag broken " , 
		     !( fabs( creal( b[ c2 + c1*NC ] ) - creal( a[ c1 + NC*c2 ] ) ) > FTOL ||
			fabs( cimag( b[ c2 + c1*NC ] ) + cimag( a[ c1 + NC*c2 ] ) ) > FTOL) ) ;
	}
      }
      //
    }
  }
  return NULL ;
}

// check if a spinor is set to the identity
static char*
identity_spinor_test( void )
{
  struct spinor S ;
  identity_spinor( &S ) ;
  size_t d1 , d2 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  // should be 1
	  if( d1 == d2 && c1 == c2 ) {
	    mu_assert( "[UNIT] error : identity_spinor broken\n" ,
		       cabs( S.D[d1][d2].C[c1][c2] - 1.0 ) < FTOL ) ;
	  // should be 0
	  } else {
	    mu_assert( "[UNIT] error : identity_spinor broken\n" ,
		       cabs( S.D[d1][d2].C[c1][c2] ) < FTOL ) ;
	  }
	}
      }
    }
  }
  return NULL ;
}

// compute B = A * Id, where Id is a color matrix, A & B are spinors
static char
*spinor_gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  spinor_gauge( &B , A , Id ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  for( i = 0 ; i < 2*NCNC*NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinor_gauge broken " , 
	       !( *a != *b ) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

// computes B = ( A )^{\dagger} Id, product over color indices
static char
*spinordag_gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ; 
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  spinordag_gauge( &B , A , Id ) ;
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      const double complex *a = (const double complex*)A.D[d1][d2].C ;
      const double complex *b = (const double complex*)B.D[d1][d2].C ;
      int c1 , c2 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  mu_assert( "[UNIT] error : spinordag_gauge broken " , 
		     !( fabs( creal( b[ c2 + c1*NC ] ) - creal( a[ c1 + NC*c2 ] ) ) > FTOL ||
			fabs( cimag( b[ c2 + c1*NC ] ) + cimag( a[ c1 + NC*c2 ] ) ) > FTOL) ) ;
	}
      }
    }
  }
  return NULL ;
}

// computes B = A * ( -I )^{\dagger}, product over color indices
static char
*spinor_gaugedag_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ; 
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : -I ;
    }
  }
  struct spinor B ;
  spinor_gaugedag( &B , A , Id ) ;
  const double complex *a = (const double complex*)A.D ;
  const double complex *b = (const double complex*)B.D ;
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor_gaugedag broken " , 
	       !( fabs( creal( b[ i ] ) + cimag( a[ i ] ) ) > FTOL ||
		  fabs( cimag( b[ i ] ) - creal( a[ i ] ) ) > FTOL) ) ;
  }
  return NULL ;
}

// sum of a timeslice test
static char
*sumprop_test( void )
{
  struct spinor SUM ;
  sumprop( &SUM , S ) ;
  const double *a = (double*)A.D ;
  const double *sum = (const double*)SUM.D ;
  int i ;
  for( i = 0 ; i < 2 * NSNS * NCNC ; i++ ) {
    mu_assert( "[UNIT] error : sum_prop broken " , 
	       !( fabs( LCU * a[i] - sum[i] ) > FTOL ) ) ;
  } 
  return NULL ;
}

// compute Id = ( Id * A ) using full spinor multiply
static char
*spinmul_atomic_left_test( void )
{
  struct spinor Id ;
  spinor_zero_site( &Id ) ; // this MUST get tested befor this function call
  // maybe we should have this as identity spinor?
  int d1d2 , c1c2 ;
  for( d1d2 = 0 ; d1d2 < NS ; d1d2++ ) {
    for( c1c2 = 0 ; c1c2 < NC ; c1c2++ ) {
      Id.D[d1d2][d1d2].C[c1c2][c1c2] = 1.0 ;
    }
  }
  spinmul_atomic_left( &Id , A ) ;
  const double *id = (const double *)Id.D ;
  const double *a = (const double *)A.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinmul_atomic_left broken " , 
	       !( fabs( *a - *id ) > FTOL ) ) ; 
    a++ , id++ ;
  }
  return NULL ;
}

// test for the whole spinor being zero
static char
*spinor_zero_test( void )
{
  spinor_zero( S ) ;
  const double *a = (const double *)S ;
  int i ;
  for( i = 0 ; i < LCU*2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor_zero broken " , 
	       !( fabs( *a ) > FTOL ) ) ; 
    a++ ;
  }
  return NULL ;
}

static char
*spintrace_test( void )
{
  spintrace( D , A.D ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    const size_t res = ( NCNC * ( NSNS - 1 ) * NS ) / 2 + i * NS ;
    mu_assert( "[UNIT] error : spinor ops spintrace broken " , 
	       !( fabs( creal( D[i] ) - res ) > FTOL 
		  ||
		  fabs( cimag( D[i] ) - res ) > FTOL 
		  ) ) ;
  }
  return NULL ;
}

// spinor tests
static char *
spinops_test( void )
{
  // initialise A
  double complex *b = (double complex*)A.D ;
  int i ;
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    *b = ( i + I * i ) , b++ ;
  }

  // set S to be A over all of LCU
  for( i = 0 ; i < LCU ; i++ ) {
    equate_spinor( &S[ i ] , A.D ) ;
  }

  // run spinor ops tests
  mu_run_test( colortrace_spinor_test ) ;
  mu_run_test( equate_minus_test ) ;
  mu_run_test( flipsign_test ) ;
  mu_run_test( gauge_test ) ;
  mu_run_test( gaugedag_test ) ;
  mu_run_test( gauge_spinordag_test ) ;
  mu_run_test( identity_spinor_test ) ;
  mu_run_test( spinor_gauge_test ) ;
  mu_run_test( spinordag_gauge_test ) ;
  mu_run_test( spinor_gaugedag_test ) ;
  mu_run_test( sumprop_test ) ;
  mu_run_test( spinmul_atomic_left_test ) ;
  mu_run_test( spinor_zero_test ) ;
  mu_run_test( spintrace_test ) ;

  return NULL ;
}

int
spinor_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;
  
  // LCU spinor storage
  S = ( struct spinor* )malloc( LCU * sizeof( struct spinor ) ) ;

  corr_malloc( (void**)&C , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&D , 16 , NCNC * sizeof( double complex ) ) ;

  // initial gamma setup and test
  char *spinres = spinops_test( ) ;

  free( D ) ;
  free( C ) ;
  free( S ) ;

  if( tests_fail == 0 ) {
    printf( "[SPINORS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , spinres ) ;
    printf( "[SPINORS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}

