/**
   @file halfspinor_tests.c
   @brief halfspinor_ops testing
 */
#include "common.h"

#include "matrix_ops.h"
#include "minunit.h"     // minimal unit testing framework
#include "halfspinor_ops.h"  // spinor operations

static struct halfspinor A ;

#define FTOL (NC*1E-14)
#define NSNS_4 (NSNS/4)

static int
halfspinor_is_identity( const struct halfspinor S )
{
  size_t i , j ;
  // test B is the identity
  for( i = 0 ; i < NC ; i++ ) {
    if( cabs( S.D[0][(NC+1)*i] - 1.0 ) > FTOL ) return FAILURE ;
    if( cabs( S.D[3][(NC+1)*i] - 1.0 ) > FTOL ) return FAILURE ;
    for( j = i+1 ; j < NC ; j++ ) {
      if( cabs( S.D[0][j+i*NC] ) > FTOL ) return FAILURE ;
      if( cabs( S.D[3][j+i*NC] ) > FTOL ) return FAILURE ;
      if( cabs( S.D[0][i+j*NC] ) > FTOL ) return FAILURE ;
      if( cabs( S.D[3][i+j*NC] ) > FTOL ) return FAILURE ;
    }
  }
  return SUCCESS ;
}

// check that (A += A) == 2A
static char*
add_halfspinor_test( void )
{
  struct halfspinor B = A ;
  add_halfspinor( &B , A ) ;
  
  const double complex *c = (const double complex*)A.D ;
  const double complex *d = (const double complex*)B.D ;
  
  size_t i ;
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops add_halfspinors broken " , 
	       !( cabs( 2*(*c) - *d ) > FTOL ) ) ;
    c++ ; d++ ;
  }
  return NULL ;
}

// check that (A*ID) == A
static char*
colormatrix_halfspinor_test( void )
{
  double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  // set M to the identity
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    C[ (NC+1)*i ] = 1.0 ;
    for( j = i+1 ; j < NC ; j++ ) {
      C[ j + i*NC ] = C[ i + j*NC ] = 0.0 ;
    }
  }

  struct halfspinor B ;
  colormatrix_halfspinor( (void*)B.D , (const void*)C , (const void*)A.D ) ;
  
  const double complex *c = (const double complex*)A.D ;
  const double complex *d = (const double complex*)B.D ;
  
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops colormatrix_halfspinor broken " , 
	       !( cabs( (*c) - *d ) > FTOL ) ) ;
    c++ ; d++ ;
  }
  return NULL ;
}

// check that (A*ID) == A
static char*
colormatrixdag_halfspinor_test( void )
{
  double complex a[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex b[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  // set M to indices
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[ i ] = 1.0 * ( 1 + I );
  }

  struct halfspinor B , D ;
  colormatrixdag_halfspinor( &B , a , A ) ;

  dagger_gauge( (void*)b , (const void*)a ) ;
  colormatrix_halfspinor( (void*)D.D , (const void*)b , (const void*)A.D ) ;
  
  const double complex *c = (const double complex*)B.D ;
  const double complex *d = (const double complex*)D.D ;
  
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops colormatrixdag_halfspinor broken " , 
	       !( cabs( (*c) - *d ) > FTOL ) ) ;
    c++ ; d++ ;
  }
  return NULL ;
}

// check that (Fmunu*A) == Fmunu*A
// by the long way
static char*
Fmunu_halfspinor_test( void )
{
  struct halfspinor B , C ;
  // create a hermitian matrix
  double complex F[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    F[ (NC+1)*i ] = 0.0 ;
    for( j = i+1 ; j < NC ; j++ ) {
      F[ j + i*NC ] = 1 + I ;
      F[ i + j*NC ] = 1 - I ;
    }
  }
  
  Fmunu_halfspinor( &B , F , A ) ;
  colormatrix_halfspinor( (void*)C.D , (const void*)F , (const void*)A.D ) ;

  const double complex *pB = (const double complex*)B.D ;
  const double complex *pC = (const double complex*)C.D ;
  
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops Fmunu_halfspinor broken " , 
    	       !( cabs( *pB - *pC ) > FTOL ) ) ;
    pB++ ; pC++ ;
  }
 
  return NULL ;
}

// check that (A += -1*A) == 0
static char*
halfspinor_Saxpy_test( void )
{
  struct halfspinor B = A ;
  halfspinor_Saxpy( &B , A , -1 ) ;

  const double complex *c = (const double complex*)B.D ;

  size_t i ;
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops halfspinor_Saxpy broken " , 
	       !( cabs( (*c) ) > FTOL ) ) ;
    c++ ;
  }
  return NULL ;
}

// check that B = A - i*A is purely real and
// equal to double the index
static char*
halfspinor_iSaxpy_test( void )
{
  struct halfspinor B = A ;
  halfspinor_iSaxpy( &B , A , -1 ) ;

  const double complex *pB = (const double complex*)B.D ;

  size_t i ;
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops halfspinor_iSaxpy broken " , 
    	       !( cabs( (*pB - 2*i ) ) > FTOL ) ) ;
    pB++ ;
  }
  return NULL ;
}

// check that B = 0 + \sigma*sigma gives the identity matrix
static char*
halfspinor_sigma_Saxpy_test( void )
{
  struct halfspinor B , S ;
  size_t i ;

  // set S to be sigma_x
  zero_halfspinor( &B ) ;
  zero_halfspinor( &S ) ;
  for( i = 0 ; i < NC ; i++ ) {
    S.D[1][(NC+1)*i] = 1 ;
    S.D[2][(NC+1)*i] = 1 ;
  }
  const uint8_t sigma_x[ NS ] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imap_x[ NS ] = { 0 , 0 , 0 , 0 } ;
  halfspinor_sigma_Saxpy( &B , S , sigma_x , imap_x ) ;

  mu_assert( "[UNIT] error : halfspinor ops halfspinor_sigma_Saxpy broken for sigmma_x" , 
	     ( halfspinor_is_identity( B ) == SUCCESS ) ) ;

  // set S to be sigma_y
  zero_halfspinor( &B ) ;
  zero_halfspinor( &S ) ;
  for( i = 0 ; i < NC ; i++ ) {
    S.D[1][(NC+1)*i] = -I ;
    S.D[2][(NC+1)*i] = +I ;
  }
  const uint8_t sigma_y[ NS ] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imap_y[ NS ] = { 3 , 3 , 1 , 1 } ;
  halfspinor_sigma_Saxpy( &B , S , sigma_y , imap_y ) ;
  mu_assert( "[UNIT] error : halfspinor ops halfspinor_sigma_Saxpy broken for sigma_y" , 
	     ( halfspinor_is_identity( B ) == SUCCESS ) ) ;

  // set S to be sigma_z
  zero_halfspinor( &B ) ;
  zero_halfspinor( &S ) ;
  for( i = 0 ; i < NC ; i++ ) {
    S.D[0][(NC+1)*i] = +1 ;
    S.D[3][(NC+1)*i] = -1 ;
  }
  const uint8_t sigma_z[ NS ] = { 0 , 1 , 2 , 3 } ;
  const uint8_t imap_z[ NS ] = { 0 , 0 , 2 , 2 } ;
  halfspinor_sigma_Saxpy( &B , S , sigma_z , imap_z ) ;
  mu_assert( "[UNIT] error : halfspinor ops halfspinor_sigma_Saxpy broken for sigma_z" , 
	     ( halfspinor_is_identity( B ) == SUCCESS ) ) ;
  
  return NULL ;
}

static char*
zero_halfspinor_test( void )
{
  struct halfspinor B ;
  zero_halfspinor( &B ) ;

  const double complex *pB = (const double complex*)B.D ;
  size_t i ;
  for( i = 0 ; i < NSNS_4*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : halfspinor ops halfspinor_iSaxpy broken " , 
    	       !( cabs( (*pB ) ) > FTOL ) ) ;
    pB++ ;
  }
  
  return NULL ;
}

// spinor tests
static char *
halfspinorops_test( void )
{
  // initialise A to be just a linear value of the indexing
  double complex *pA = (double complex*)A.D ; 
  size_t i ;
  for( i = 0 ; i < (NSNS_4)*NCNC ; i++ ) {
    *pA = i*( 1 + I ) ; pA++ ;
  }

  mu_run_test( add_halfspinor_test ) ;
  mu_run_test( colormatrix_halfspinor_test ) ;
  // relies on colormatrix_halfspinor to work
  mu_run_test( Fmunu_halfspinor_test ) ;
  mu_run_test( colormatrixdag_halfspinor_test ) ;
  mu_run_test( halfspinor_Saxpy_test ) ;
  mu_run_test( halfspinor_iSaxpy_test ) ;
  mu_run_test( zero_halfspinor_test ) ;
  // relies on zero_halfspinor
  mu_run_test( halfspinor_sigma_Saxpy_test ) ;
  
  return NULL ;
}

int
halfspinor_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;
  
  // initial gamma setup and test
  char *spinres = halfspinorops_test( ) ;

  if( tests_fail == 0 ) {
    fprintf( stdout , "[HALFSPINORS UNIT] all %d tests passed\n\n" ,
	     tests_run ) ;
    return SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , spinres ) ;
    fprintf( stderr , "[HALFSPINORS UNIT] %d out of %d tests failed\n\n" , 
	     tests_fail , tests_run ) ;
    return FAILURE ;
  }
}

#undef FTOL
