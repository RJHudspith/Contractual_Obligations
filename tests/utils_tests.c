/**
   @file utils_tests.c
   @brief test of the codes in the UTILS directory
   @warning must be called after spinor tests to ensure identity spinor checked
 */
#include "common.h"

#include "basis_conversions.h"  // conversion between chiral and *
#include "corr_malloc.h"        // test we can allocate stuff
#include "corr_sort.h"          // test our sort code
#include "crc32c.h"             // test our crc
#include "GLU_bswap.h"          // test the byte swap
#include "gramschmidt.h"        // test our gramschmidt code
#include "minunit.h"            // unit test framework
#include "spinor_ops.h"         // identity_spinor()

#define FTOL (1E-14*NC)

// test that rotation of identity is still identity
static char *
basis_conversion_test( void )
{
  // basis rotation of identity spinor should be identity spinor
  struct spinor Id , Rt ;
  identity_spinor( &Id ) ;
  Rt = Id ;
  chiral_to_nrel( &Rt ) ;
  size_t d1 , d2 , c1, c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  mu_assert( "[UNIT] error : basis conversion failed\n" ,
		     cabs( Id.D[d1][d2].C[c1][c2] 
			   - Rt.D[d1][d2].C[c1][c2] ) < 1E-14 ) ;
	}
      }
    }	
  }
  return NULL ;
}

// little checking code
static uintptr_t is_aligned( const void *a , 
			     const size_t size ) 
{
  return (uintptr_t)a%size == 0 ;
}

// check our malloc
static char *
corr_malloc_test( void )
{
  int *a = NULL ;
  const int errorcode = corr_malloc( (void*)&a , 16 , 1 ) ;
  // check we can actually allocate something
  mu_assert( "[UNIT] error : cannot allocate memory\n" ,
	     errorcode == 0 ) ;
  // check it is aligned
#ifdef HAVE_IMMINTRIN_H
  mu_assert( "[UNIT] error : memory is not aligned\n" ,
	     is_aligned( a , 16 ) ) ;
#endif
  free( a ) ;
  return NULL ;
}

// test our crc32c against known input
static char *
crc32c_test( void )
{
  // test against zero
  uint32_t a[1] = { 0x00000000 } , ca = 0 , cb = 0 ;
  DML_checksum_accum_crc32c( &ca , &cb , 0 , a , sizeof( uint32_t ) ) ;
  mu_assert( "[UNIT] error : crc32c broken\n" ,
	     ( ca == 0x00000000 && cb == 0x00000000 ) ) ;

  // test against max int
  a[0] = 0xFFFFFFFF , ca = 0 , cb = 0 ;
  DML_checksum_accum_crc32c( &ca , &cb , 32 , a , sizeof( uint32_t ) ) ;
  mu_assert( "[UNIT] error : crc32c broken\n" ,
	     ( ca == 0xBCC5A1C5 && cb == 0x6F316871 ) ) ;

  return NULL ;
}

// union to test byte swaps
union uint32_dbl{
  uint32_t u[ 2 ] ;
  double val ;
} ;

// test our byte swaps
static char *
GLU_bswap_test( void )
{
  uint16_t a = 0x00FF ;
  bswap_16( 1 , &a ) ;
  mu_assert( "[UNIT] error : bswap_16 broken \n" ,
	     ( a == 0xFF00 ) ) ;
  uint32_t b = 0x0000FFFF ;
  bswap_32( 1 , &b ) ;
  mu_assert( "[UNIT] error : bswap_32 broken \n" ,
	     ( b == 0xFFFF0000 ) ) ;
  union uint32_dbl c = { .u[0] = 0x00000000 ,
			 .u[1] = 0x0000FFFF } ;
  double cc = c.val ;
  bswap_64( 1 , &cc ) ;
  c.u[0] = 0xFFFF0000 ; c.u[1] = 0x00000000 ;
  mu_assert( "[UNIT] error : bswap_64 broken \n" ,
	     ( c.val == cc ) ) ;
  return NULL ;
}

// test our gram-schmidt routine
static char *
gramschmidt_test( void )
{
  double complex a[ NCNC ] ;
  size_t i , j , k ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] = i * ( 1 + I ) ;
  }
  reunit2( a ) ;

  // test a * a^\dagger = Id
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      double complex sum = 0.0 ;
      for( k = 0 ; k < NC ; k++ ) {
	sum += a[ k + j*NC ] * conj( a[ k + i*NC ] ) ;
      }
      if( i == j ) {
	mu_assert( "[UNIT] error : gramschmidt broken \n" ,
		   ( cabs( 1.0 - sum ) < FTOL ) ) ;
      } else {
	mu_assert( "[UNIT] error : gramschmidt broken \n" ,
		   ( cabs( sum ) < FTOL ) ) ;
      }
    }
  }

  return NULL ;
}

// test that a descending list of ints becomes an ascending one
static char *
merge_sort_test( void )
{
  const size_t N = 1E2 ;
  int *a = malloc( N * sizeof( int ) ) ;
  int *b = malloc( N * sizeof( int ) ) ;
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    a[i] = b[i] = N-i-1 ;
  }
  merge_sort( a , b , sizeof(int) , sizeof(int) ,
	      N , lt_int ) ;
  for( i = 0 ; i < N ; i++ ) {
    mu_assert( "[UTILS] error : merge_sort broken" , 
	       ( a[i] == (int)i && b[i] == (int)i ) ) ;
  }
  return NULL ;
}

// run the utils tests
static char *
utils_test( void )
{
  // check the basis conversion
  mu_run_test( basis_conversion_test ) ;
  mu_run_test( crc32c_test ) ;
  mu_run_test( corr_malloc_test ) ;
  mu_run_test( GLU_bswap_test ) ;
  mu_run_test( gramschmidt_test ) ;
  mu_run_test( merge_sort_test ) ;
  return NULL ;
}

// full tests
int
utils_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  char *utils_res = utils_test( ) ;

  if( tests_fail == 0 ) {
    fprintf( stdout , "[UTILS UNIT] all %d tests passed\n\n" ,
	     tests_run ) ;
    return SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , utils_res ) ;
    fprintf( stderr , "[UTILS UNIT] %d out of %d tests failed\n\n" , 
	     tests_fail , tests_run ) ;
    return FAILURE ;
  }

}
