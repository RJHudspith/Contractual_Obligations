/**
   @file matops_tests.c
   @brief matrix operations tests
*/

#include "common.h"

#include "minunit.h"     // minimal unit testing framework
#include "matrix_ops.h"  // matrix operations

#ifdef HAVE_IMMINTRIN_H
typedef __m128d dcomplex ;
#else
typedef double complex dcomplex ;
#endif

static double complex *U  ; // is a color matrix
static double complex *V  ; // is a color matrix
static double complex *Id ; // is a color matrix
static dcomplex *res ;

#define FTOL ( NC * 1.E-14 ) 

// needed for casting a possible __m128d to double complex
static inline double complex
dcast( dcomplex a )
{
#ifdef HAVE_IMMINTRIN_H
  double complex s ;
  _mm_store_pd( (void*)&s , a ) ;
  return s ;
#else
  return a ;
#endif
}

// tests U + U = 2 U
static char *
add_mat_test( void )
{
  colormatrix_equiv( V , U ) ;
  add_mat( (dcomplex*)V , (const dcomplex*)U ) ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    mu_assert( "[MATOPS UNIT] error : add_mat broken" , 
	       !( fabs( creal( V[i] ) - 2.0 * creal( U[i] ) ) ||
		  fabs( cimag( V[i] ) - 2.0 * cimag( U[i] ) ) ) ) ;
  }
  return NULL ;
}

// tests equiv U - V = 0
static char *
colormatrix_equiv_test( void )
{
  double complex V[ NCNC ] ;
  colormatrix_equiv( V , U ) ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    mu_assert( "[MATOPS UNIT] error : colormatrix_equiv broken" , 
	       !( fabs( creal( V[i] ) - creal( U[i] ) ) > FTOL ||
		  fabs( cimag( V[i] ) - cimag( U[i] ) ) > FTOL ) ) ;
  }
  return NULL ;
}

// tests equiv U - V = 0
static char *
colormatrix_equiv_f2d_test( void )
{
  float complex v[ NCNC ] ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    v[ i ] = (float complex)U[ i ] ;
  }
  double complex W[ NCNC ] ;
  colormatrix_equiv_f2d( W , v ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    mu_assert( "[MATOPS UNIT] error : colormatrix_equiv broken" , 
	       !( fabs( creal( W[i] ) - creal( U[i] ) ) > 1E-6 ||
		  fabs( cimag( W[i] ) - cimag( U[i] ) ) > 1E-6 ) ) ;
  }
  return NULL ;
}

// trace of a matrix
static char *
colortrace_test( void )
{
  double complex tr = dcast( colortrace( (const dcomplex*)U ) ) ;
  // closed form solution
  const double sol = 0.5 * ( NC * ( NC - 1 ) ) * ( NC + 1 ) ;
  mu_assert( "[MATOPS UNIT] error : colortrace broken" , 
	     !( fabs( creal( tr ) - sol ) > FTOL ||
		fabs( cimag( tr ) - sol ) > FTOL ) ) ;
  return NULL ;
}

// trace of U^2
static char *
colortrace_prod_test( void )
{
  double complex tr = dcast( colortrace_prod( (const dcomplex*)U , 
					      (const dcomplex*)U ) ) ;
  multab( res , (const dcomplex*)U , (const dcomplex*)U ) ;
  double complex tr2 = dcast( colortrace( res ) ) ;
  // need some closed form solution here
  mu_assert( "[MATOPS UNIT] error : colortrace_prod broken" , 
	     !( fabs( creal( tr ) - creal( tr2 ) ) > FTOL ||
		fabs( cimag( tr ) - cimag( tr2 ) ) > FTOL ) ) ;
  return NULL ;
}

// test the conjugate transpose
static char *
dagger_test( void )
{
  dagger_gauge( res , (dcomplex*)U ) ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      const double complex r = dcast( res[ j + i*NC ] ) ;
      const double complex t = U[ i + j*NC ] ;
      mu_assert( "[MATOPS UNIT] error : dagger broken" , 
		 !( fabs( creal( r ) - creal( t ) ) > FTOL ||
		    fabs( cimag( r ) + cimag( t ) ) > FTOL ) ) ;
    }
  }
  return NULL ;
}

// matrix multiply - Multiply by the identity as a first test
static char *
multab_test( void )
{
  int i ;
  multab( res , (const dcomplex*)Id , (const dcomplex*)U ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    const double complex r = dcast( res[ i ] ) ;
    mu_assert( "[MATOPS UNIT] error : multab broken" , 
	       !( fabs( creal( r ) - creal( U[i] ) ) > FTOL ||
		  fabs( cimag( r ) - cimag( U[i] ) ) > FTOL ) ) ;
  }
  return NULL ;
}

// multiply by the identity and dagger U
static char *
multab_dag_test( void )
{
  int i , j ;
  multab_dag( res , (const dcomplex*)Id , (const dcomplex*)U ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      const double complex r = dcast( res[ j + i*NC ] ) ;
      mu_assert( "[MATOPS UNIT] error : multab_dag broken" , 
		 !( fabs( creal( r ) - creal( U[ i + j*NC ] ) ) > FTOL ||
		    fabs( cimag( r ) + cimag( U[ i + j*NC ] ) ) > FTOL ) ) ;
    }
  }
  return NULL ;
}

// dagger U and multiply by the identity
static char *
multabdag_test( void )
{
  int i , j ;
  multabdag( res , (const dcomplex*)U , (const dcomplex*)Id ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      const double complex r = dcast( res[ j + i*NC ] ) ;
      mu_assert( "[MATOPS UNIT] error : multabdag broken" , 
		 !( fabs( creal( r ) - creal( U[ i + j*NC ] ) ) > FTOL ||
		    fabs( cimag( r ) + cimag( U[ i + j*NC ] ) ) > FTOL ) ) ;
    }
  }
  return NULL ;
}

// dagger U and multiply by the identity
static char *
multabdagdag_test( void )
{
  int i , j ;
  multab_dagdag( res , (const dcomplex*)U , (const dcomplex*)Id ) ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      const double complex r = dcast( res[ j + i*NC ] ) ;
      mu_assert( "[MATOPS UNIT] error : multabdagdag broken" , 
		 !( fabs( creal( r ) - creal( U[ i + j*NC ] ) ) > FTOL ||
		    fabs( cimag( r ) + cimag( U[ i + j*NC ] ) ) > FTOL ) ) ;
    }
  }
  return NULL ;
}

// spinor tests
static char *
matops_test( void )
{
  // gauge matrix
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    U[ i ] = ( i + I * i ) ; 
  }

  // identity matrix
  int j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( i != j ) ? 0.0 : 1.0 ;
    }
  }
  
  // tests
  mu_run_test( colormatrix_equiv_test ) ;
  mu_run_test( colormatrix_equiv_f2d_test ) ;
  mu_run_test( add_mat_test ) ;
  mu_run_test( colortrace_test ) ;
  mu_run_test( dagger_test ) ;
  mu_run_test( multab_test ) ;
  mu_run_test( multab_dag_test ) ;
  mu_run_test( multabdag_test ) ;
  mu_run_test( multabdagdag_test ) ;
  mu_run_test( colortrace_prod_test ) ;

  return NULL ;
}

// driver for the tests
int
matops_test_driver( void )
{
  // temporary storage
  corr_malloc( (void**)&U   , 16 , NCNC * sizeof( double complex ) ) ;
  corr_malloc( (void**)&V   , 16 , NCNC * sizeof( double complex ) ) ;
  corr_malloc( (void**)&res , 16 , NCNC * sizeof( double complex ) ) ;
  corr_malloc( (void**)&Id  , 16 , NCNC * sizeof( double complex ) ) ;

  // initialise the test counters
  tests_run = tests_fail = 0 ;

  // matrix operations test
  char *matres = matops_test( ) ;

  // and free the memory
  free( U ) ;
  free( V ) ;
  free( res ) ;
  free( Id ) ;

  if( tests_fail == 0 ) {
    fprintf( stdout , "[MATOPS UNIT] all %d tests passed\n\n" ,
	     tests_run ) ;
    return SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , matres ) ;
    fprintf( stderr , "[MATOPS UNIT] %d out of %d tests failed\n\n" , 
	     tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
