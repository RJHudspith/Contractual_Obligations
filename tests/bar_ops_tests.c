/**
   @file bar_ops_tests.c
   @brief baryon operations tests
 */
#include "common.h"

#include "bar_ops.h"
#include "minunit.h"

// our tolerance
#define FLTOL (NC*1.E-14)

// some temporary space
static struct spinor S1 , S2 ;

// check the baryon contractor
static char *
baryon_contract_test( void )
{
  // closed-form solution for each trace
  const size_t res = NC * ( NC - 1 ) * NC * ( 2*NC - 1 ) / 3 ;
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      const double complex r = baryon_contract( S1 , S2 ,
						d1 , d2 , 
						d1 , d2 ) ;
      // should be purely imaginary
      mu_assert( "[UNIT] error : bar_ops baryon_contract broken\n" ,
		 fabs( creal( r ) + cimag( r ) - res ) < FLTOL ) ;
    }
  }
  return NULL ;
}

// check the cross color product: should be zero with our
// color content as det(c) == 0
static char *
cross_color_test( void )
{
  cross_color_trace( &S1 , S2 ) ;
  const double complex *a = (const double complex*)S1.D ;
  size_t d ;
  for( d = 0 ; d < NSNS*NCNC ; d++ ) {
    mu_assert( "[UNIT] error : bar_ops cross_color_trace broken\n" ,
	       cabs( a[ d ] ) < FLTOL ) ;
  }
  return NULL ;
}

// baryon operations tests
static char *
bar_ops_test( void )
{
  // initialise spinors
  double complex *a = (double complex*)S1.D ;
  double complex *b = (double complex*)S2.D ;
  double complex c[ NCNC ] ; // special color matrix det(c) = 0!
  size_t i , j ;
  // init the color matrix
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      c[ j + i * NC ] = ( i + I * i ) ;
    }
  }
  // set up the spinor
  for( i = 0 ; i < NSNS ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      *b = c[ j ] ;
      *a = *b ;
      a++ , b++ ;
    }
  }

  // check the baryon contraction
  mu_run_test( baryon_contract_test ) ;

  // check the cross product is 0
  mu_run_test( cross_color_test ) ;

  return NULL ;
}

// runs the whole #!
int
bar_ops_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  // initial gamma setup and test
  char *baropsres = bar_ops_test( ) ;

  if( tests_fail == 0 ) {
    printf( "[BAR_OPS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , baropsres ) ;
    printf( "[BAR_OPS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
