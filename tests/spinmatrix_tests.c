/**
   @file spinor_tests.c
   @brief spinor_ops testing
 */
#include "common.h"

#include "gammas.h"
#include "minunit.h"     // minimal unit testing framework
#include "spinmatrix_ops.h"  // spinor operations

static struct gamma *GAMMA = NULL ; // gamma matrix technology

static double complex *D = NULL ;
static double complex *C = NULL ;

#define FTOL ( NC * 1.E-14 ) 

// add spinmatrix test
static char*
atomic_add_spinmatrices_test( void )
{
  // set C to be the identity
  identity_spinmatrix( C ) ;
  atomic_add_spinmatrices( C , D ) ;
  // check that C = 1 + D
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      mu_assert( "[UNIT] error : spinmatrix ops atomic_add_spinmatrices broken " , 
		 !( i == j ? \
		   cabs( C[ j + i * NS ] - 1.0 - D[ j + i * NS ] ) > FTOL : \
		   cabs( C[ j + i * NS ] - D[ j + i * NS ] ) > FTOL )
		 ) ;
    }
  }
  return NULL ;
}

// gamma spinmatrix test, multiply by identity?
static char*
gamma_spinmatrix_test( void )
{
  // set c to d and multiply c by identity
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  gamma_spinmatrix( C , GAMMA[ IDENTITY ] ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops gamma_spinmatrix broken " , 
	       !( cabs( C[ i ] - D[ i ] ) > FTOL ) ) ;
  }
  return NULL ;
}

// gammaspin_trace test
static char*
gammaspinmatrix_trace_test( void )
{
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  // direct trace of the product
  const double complex tr1 = gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) ;
  gamma_spinmatrix( C , GAMMA[ GAMMA_3 ] ) ;
  const double complex tr2 = spinmatrix_trace( C ) ; // trace( mul() )
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  spinmatrix_gamma( C , GAMMA[ GAMMA_3 ] ) ;
  const double complex tr3 = spinmatrix_trace( C ) ; // cyclicity
  mu_assert( "[UNIT] error : spinmatrix ops gammaspinmatrix_trace broken " , 
	       !( cabs( tr1 - tr2 ) > FTOL ) ) ;
  mu_assert( "[UNIT] error : spinmatrix ops gammaspinmatrix_trace broken " , 
	       !( cabs( tr1 - tr3 ) > FTOL ) ) ;
  return NULL ;
}

// trace should be NS
static char*
identity_spinmatrix_test( void )
{
  identity_spinmatrix( C ) ;
  const double complex tr = spinmatrix_trace( C ) ;
  mu_assert( "[UNIT] error : spinmatrix ops identity_spinmatrix broken " , 
	     !( cabs( tr - NS ) > FTOL ) ) ;
  return NULL ;
}

// spinmatrix gamma test, multiply by identity?
static char*
spinmatrix_gamma_test( void )
{
  // set c to d and multiply c by identity
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  spinmatrix_gamma( C , GAMMA[ IDENTITY ] ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_gamma broken " , 
	       !( cabs( C[ i ] - D[ i ] ) > FTOL ) ) ;
  }
  return NULL ;
}

// multiply by zero test
static char*
spinmatrix_mulconst_test( void )
{
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  spinmatrix_mulconst( C , 0.0 ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_mulconst broken " , 
	       !( cabs( C[ i ] ) > FTOL ) ) ;
  }
  return NULL ;
}

// multiply
static char*
spinmatrix_multiply_test( void )
{
  identity_spinmatrix( C ) ;
  double complex A[ NSNS ] ;
  spinmatrix_multiply( A , C , D ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_multiply broken " , 
	       !( cabs( A[i] - D[i] ) > FTOL ) ) ;
  }
  return NULL ;
}

// trace test uses an identity
static char*
spinmatrix_trace_test( void )
{
  const double complex tr = spinmatrix_trace( D ) ;
  const int sol = NS * ( NSNS - 1 ) / 2 ;
  mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_trace broken " , 
	     !( fabs( creal( tr ) - sol ) > FTOL || 
		fabs( cimag( tr ) - sol ) > FTOL ) ) ;
  return NULL ;
}

// test if we can zero an array
static char*
zero_spinmatrix_test( void )
{
  zero_spinmatrix( C ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops zero_spinmatrix broken " , 
	       !( cabs( C[i] ) > FTOL ) ) ;
  }
  return NULL ;
}

// spinor tests
static char *
spinmatrices_test( void )
{
  // initialise D
  double complex *b = (double complex*)D ;
  int i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *b = ( i + I * i ) , b++ ;
  }

  // run spinmatrix ops tests
  mu_run_test( gamma_spinmatrix_test ) ;
  mu_run_test( identity_spinmatrix_test ) ;
  mu_run_test( atomic_add_spinmatrices_test ) ;
  mu_run_test( spinmatrix_gamma_test ) ;
  mu_run_test( spinmatrix_mulconst_test ) ;
  mu_run_test( spinmatrix_multiply_test ) ;
  mu_run_test( spinmatrix_trace_test ) ;
  mu_run_test( zero_spinmatrix_test ) ;

  // relies on spinmatrix multiply and spinmatrix trace
  mu_run_test( gammaspinmatrix_trace_test ) ;

  return NULL ;
}

// runs the whole #!
int
spinmatrix_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;
  
  // precompute the gamma basis
  GAMMA = malloc( NSNS * sizeof( struct gamma ) ) ;

  corr_malloc( (void**)&C , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&D , 16 , NSNS * sizeof( double complex ) ) ;

  make_gammas( GAMMA , CHIRAL ) ;

  // initial gamma setup and test
  char *spinres = spinmatrices_test( ) ;

  // memfree
  free( C ) ;
  free( D ) ;
  free( GAMMA ) ;

  if( tests_fail == 0 ) {
    printf( "[SPINMATRIX UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , spinres ) ;
    printf( "[SPINMATRIX UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}

