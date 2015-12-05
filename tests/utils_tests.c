/**
   @file utils_tests.c
   @brief test of the codes in the UTILS directory
   @warning must be called after spinor tests to ensure identity spinor checked
 */
#include "common.h"

#include "basis_conversions.h"  // conversion between chiral and *
#include "minunit.h"            // unit test framework
#include "spinor_ops.h"         // identity_spinor()

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

// run the utils tests
static char *
utils_test( void )
{
  // check the basis conversion
  mu_run_test( basis_conversion_test ) ;

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
    printf( "[UTILS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , utils_res ) ;
    printf( "[UTILS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }

}
