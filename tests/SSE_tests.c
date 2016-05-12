/**
   @file SSE_tests.c
   @brief check our SSE'd complex multiplies
 */
#include "common.h"
#include "minunit.h"
#include "SSE2_OPS.h"

#define FL_TOL (1E-14)

#ifdef HAVE_EMMINTRIN_H

// tests I*a SSE2/3 variant
static char *
SSE2_CONJ_test( void )
{
  double complex a = 1 + I , b ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  _mm_store_pd( (void*)&b , SSE2_CONJ( A ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_miMUL broken" , 
	     ( cabs( conj(a) - b ) < FL_TOL ) ) ;
  return NULL ;
}

// test we can flip the sign of the real and imaginary parts
static char *
SSE_FLIP_test( void )
{
  double complex a = 1 + I , b = 0.0 ;
  __m128d A = _mm_set_pd( 1.0 , 1.0 ) ;
  // cast back to double complex
  _mm_store_pd( (void*)&b , SSE_FLIP( A ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_FLIP broken" , 
	     ( cabs( a + b ) < FL_TOL ) ) ;
  return NULL ;
}

// tests I*a SSE2/3 variant
static char *
SSE2_iMUL_test( void )
{
  double complex a = 1 + I , b ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  _mm_store_pd( (void*)&b , SSE2_iMUL( A )  ) ;
  mu_assert( "[SSE_OPS] error : SSE_iMUL broken" , 
	     ( cabs( I*a - b ) < FL_TOL ) ) ;
  return NULL ;
}

// tests I*a SSE2/3 variant
static char *
SSE2_miMUL_test( void )
{
  double complex a = 1 + I , b ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  _mm_store_pd( (void*)&b , SSE2_miMUL( A )  ) ;
  mu_assert( "[SSE_OPS] error : SSE_miMUL broken" , 
	     ( cabs( -I*a - b ) < FL_TOL ) ) ;
  return NULL ;
}

// tests a*b SSE2/3 variant
static char *
SSE_MUL_test( void )
{
  double complex a = 1 + I , b = 1 + 2 * I , c ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  __m128d B = _mm_setr_pd( creal( b ) , cimag( b ) ) ;
  // cast back to double complex
  _mm_store_pd( (void*)&c , SSE2_MUL( A , B ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_MUL broken" , 
	     ( cabs( a*b - c ) < FL_TOL ) ) ;
  return NULL ;
}

// tests conj(a)*b SSE2/3 variant
static char *
SSE_MULCONJ_test( void )
{
  double complex a = 1 + I , b = -1 + 2 * I , c ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  __m128d B = _mm_setr_pd( creal( b ) , cimag( b ) ) ;
  // cast back to double complex
  _mm_store_pd( (void*)&c , SSE2_MULCONJ( A , B ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_MULCONJ broken" , 
	     ( cabs( conj(a)*b - c ) < FL_TOL ) ) ;
  return NULL ;
}

// tests a*conj(b) SSE2/3 variant
static char *
SSE_MUL_CONJ_test( void )
{
  double complex a = 1 + I , b = -1 + 2 * I , c ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  __m128d B = _mm_setr_pd( creal( b ) , cimag( b ) ) ;
  // cast back to double complex
  _mm_store_pd( (void*)&c , SSE2_MUL_CONJ( A , B ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_MUL_CONJ broken" , 
	     ( cabs( a*conj(b) - c ) < FL_TOL ) ) ;
  return NULL ;
}

// tests conj(a)*conj(b) SSE2/3 variant
static char *
SSE_MUL_CONJCONJ_test( void )
{
  double complex a = 1 + I , b = -1 + 2 * I , c ;
  __m128d A = _mm_setr_pd( creal( a ) , cimag( a ) ) ;
  __m128d B = _mm_setr_pd( creal( b ) , cimag( b ) ) ;
  // cast back to double complex
  _mm_store_pd( (void*)&c , SSE2_MUL_CONJCONJ( A , B ) ) ;
  mu_assert( "[SSE_OPS] error : SSE_MUL_CONJCONJ broken" , 
	     ( cabs( conj(a)*conj(b) - c ) < FL_TOL ) ) ;
  return NULL ;
}

#endif

// run the utils tests
static char *
SSE_OPS_test( void )
{
#ifdef HAVE_EMMINTRIN_H
  // check the basis conversion
  mu_run_test( SSE2_CONJ_test ) ;
  mu_run_test( SSE_FLIP_test ) ;
  mu_run_test( SSE2_iMUL_test ) ;
  mu_run_test( SSE2_miMUL_test ) ;
  mu_run_test( SSE_MUL_test ) ;
  mu_run_test( SSE_MULCONJ_test ) ;
  mu_run_test( SSE_MUL_CONJ_test ) ;
  mu_run_test( SSE_MUL_CONJCONJ_test ) ;
#endif
  return NULL ;
}

// full tests
int
SSE_OPS_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  char *SSE_OPS_res = SSE_OPS_test( ) ;

  if( tests_fail == 0 ) {
    fprintf( stdout , "[SSE_OPS UNIT] all %d tests passed\n\n" ,
	     tests_run ) ;
    return SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , SSE_OPS_res ) ;
    fprintf( stderr , "[SSE_OPS UNIT] %d out of %d tests failed\n\n" , 
	     tests_fail , tests_run ) ;
    return FAILURE ;
  }
}

#undef FL_TOL
