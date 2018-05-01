/**
   @file spinmatrix_tests.c
   @brief spinmatrix operations testing
 */
#include "common.h"

#include "gammas.h"
#include "minunit.h"         // minimal unit testing framework
#include "spinmatrix_ops.h"  // spinmatrix operations
#include "spinor_ops.h"      // identity_spinor()

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

// multiply by 1 should be the sum of our gamma matrices
static char*
compute_pslash_test( void ) 
{
  // set a spinmatrix to the identity
  zero_spinmatrix( D ) ;
  double p[ ND ] ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    p[ mu ] = 1.0 ;
    identity_spinmatrix( C ) ;
    gamma_spinmatrix( C , GAMMA[ mu ] ) ;
    atomic_add_spinmatrices( D , C ) ;
  }
  compute_pslash( C , GAMMA , p ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinmatrix ops compute_pslash broken " , 
	       !( cabs( C[i] - D[i] ) > FTOL ) ) ;
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
  const double complex tr1 = gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) ;
  gamma_spinmatrix( C , GAMMA[ GAMMA_T ] ) ;
  const double complex tr2 = spinmatrix_trace( C ) ; // trace( mul() )
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  spinmatrix_gamma( C , GAMMA[ GAMMA_T ] ) ;
  const double complex tr3 = spinmatrix_trace( C ) ; // cyclicity
  mu_assert( "[UNIT] error : spinmatrix ops gammaspinmatrix_trace broken " , 
	       !( cabs( tr1 - tr2 ) > FTOL ) ) ;
  mu_assert( "[UNIT] error : spinmatrix ops gammaspinmatrix_trace broken " , 
	       !( cabs( tr1 - tr3 ) > FTOL ) ) ;
  return NULL ;
}

// test our routine for peeking spin indices
static char*
get_spinmatrix_test( void )
{
  struct spinor Id ;
  double complex s[ NSNS ] ;
  identity_spinor( &Id ) ;
  size_t c1 , c2 , d1d2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      get_spinmatrix( s , Id , c1 , c2 ) ;
      // check it is the identity
      for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
	if( d1d2%(NS+1)==0 && c1==c2 ) {
	  mu_assert( "[UNIT] error : spinmatrix ops get_spinmatrix broken " , 
		     !( cabs( s[d1d2] - 1.0 ) > FTOL ) ) ;
	} else {
	  mu_assert( "[UNIT] error : spinmatrix ops get_spinmatrix broken " , 
		     !( cabs( s[d1d2] ) > FTOL ) ) ;
	}
	//
      }
      //
    }
  }
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

// square a matrix and take its trace
static char*
spinmatrix_multiply_test( void )
{
  double complex A[ NSNS ] __attribute__((aligned(ALIGNMENT)));

  // compute C = D^{T}
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      C[ d2 + d1*NS ] = D[ d1 + d2*NS ] ;
    }
  }

  // multiply D with its transpose
  spinmatrix_multiply( A , D , C ) ;

  // do a slow, by hand matrix multiply to check
  register double complex sum ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      size_t d3 ;
      sum = 0 ;
      for( d3 = 0 ; d3 < NS ; d3++ ) {
	sum += D[ d3 + d2 * NS ] * D[ d3 + d1 * NS ] ;
      }
      mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_multiply broken " , 
		 !( cabs( sum - A[d2 + d1*NS] ) > PREC_TOL ) ) ;
    }
  }
 
  return NULL ;
}

// test against transpose and multiply
static char*
spinmatrix_multiply_T_test( void )
{
  spinmatrix_multiply_T( C , D , D ) ;

  double complex A[ NSNS ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NSNS ] __attribute__((aligned(ALIGNMENT))) ;

  memcpy( A , D , NSNS * sizeof( double complex ) ) ;
  transpose_spinmatrix( A ) ;

  spinmatrix_multiply( B , D , A ) ;

  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
      mu_assert( "[UNIT] error : spinmatrix ops spinmatrix_multiply_T broken " , 
		 !( cabs( B[i] - C[i] ) > PREC_TOL ) ) ;
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

// trace of the product of two matrices
static char*
trace_prod_spinmatrices_test( void )
{
  const double complex tr1 = trace_prod_spinmatrices( D , D ) ;
  spinmatrix_multiply( C , D , D ) ;
  const double complex tr2 = spinmatrix_trace( C ) ;
  mu_assert( "[UNIT] error : spinmatric ops trace_prod_spinmatrices"
	     " broken" , !( cabs( tr1 - tr2 ) > FTOL ) ) ;
  return NULL ;
}

// test we can do a transpose right
static char*
transpose_spinmatrix_test( void )
{
  memcpy( C , D , NSNS * sizeof( double complex ) ) ;
  transpose_spinmatrix( C ) ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      mu_assert( "[UNIT] error : spinmatric ops transpose_spinmatrix"
		 " broken" , !( cabs( C[j+i*NS] - D[i+j*NS] ) > FTOL ) ) ;
    }
  }
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

// spinmatrix tests
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
  mu_run_test( get_spinmatrix_test ) ;
  mu_run_test( identity_spinmatrix_test ) ;
  mu_run_test( atomic_add_spinmatrices_test ) ;
  mu_run_test( spinmatrix_gamma_test ) ;
  mu_run_test( spinmatrix_mulconst_test ) ;
  mu_run_test( spinmatrix_multiply_test ) ;
  mu_run_test( spinmatrix_trace_test ) ;
  mu_run_test( transpose_spinmatrix_test ) ;
  mu_run_test( zero_spinmatrix_test ) ;

  // relies on spinmatrix multiply and spinmatrix trace
  mu_run_test( trace_prod_spinmatrices_test ) ;
  mu_run_test( gammaspinmatrix_trace_test ) ;
  mu_run_test( compute_pslash_test ) ;

  // relies on transpose and spinmatrix multiply
  mu_run_test( spinmatrix_multiply_T_test ) ;
  
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

  corr_malloc( (void**)&C , ALIGNMENT , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&D , ALIGNMENT , NSNS * sizeof( double complex ) ) ;

  make_gammas( GAMMA , CHIRAL ) ;

  // initial gamma setup and test
  char *spinres = spinmatrices_test( ) ;

  // memfree
  free( C ) ;
  free( D ) ;
  free( GAMMA ) ;

  if( tests_fail == 0 ) {
    fprintf( stdout , "[SPINMATRIX UNIT] all %d tests passed\n\n" ,
	     tests_run ) ;
    return SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , spinres ) ;
    fprintf( stderr , "[SPINMATRIX UNIT] %d out of %d tests failed\n\n" , 
	     tests_fail , tests_run ) ;
    return FAILURE ;
  }
}

