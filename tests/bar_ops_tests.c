/**
   @file bar_ops_tests.c
   @brief baryon operations and projections tests
 */
#include "common.h"

#include "bar_contractions.h"
#include "bar_ops.h"
#include "gammas.h"
#include "minunit.h"
#include "matrix_ops.h"
#include "spinor_ops.h" // spinor_identity
#include "spinmatrix_ops.h"

// our tolerance
#define FLTOL (NC*1.E-14)

// some temporary space
static struct spinor S1 , S2 ;

static char*
baryon_contract_site_test( void )
{
  // 3 temporary spinors
  struct spinor a , b , c ;
  spinor_zero_site( &a ) ;
  spinor_zero_site( &b ) ;
  spinor_zero_site( &c ) ;

  // lower triangular matrix all ones, determinant is 1 so
  // cross color trace is simply identity matrix multiplied 
  // by 2 * NC
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    size_t i , j ;
    for( i = 0 ; i < NC ; i++ ) {
      for( j = 0 ; j <= i ; j++ ) {
	a.D[d][d].C[i][j] =			\
	  b.D[d][d].C[i][j] =			\
	  c.D[d][d].C[i][j] = 1.0 ;
      }
    }
  }

  // allocate and set gammas
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  make_gammas( GAMMAS , CHIRAL ) ;

  // allocate and set the terms
  double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
  term[0] = malloc( NSNS * sizeof( double complex ) ) ;
  term[1] = malloc( NSNS * sizeof( double complex ) ) ;

  size_t GSGK ;
  for( GSGK = 0 ; GSGK < (B_CHANNELS*B_CHANNELS) ; GSGK++ ) {

    // zero the terms
    zero_spinmatrix( term[0] ) ;
    zero_spinmatrix( term[1] ) ;

    // compute the gammas we will be looking at
    const size_t GSRC = GSGK / B_CHANNELS ;
    const size_t GSNK = GSGK % B_CHANNELS ;
    const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
    const struct gamma Cgnu = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
    const struct gamma CgnuT = CGmuT( Cgnu , GAMMAS ) ;
    baryon_contract_site( term , a ,  b , c , Cgmu , CgnuT ) ;

    // compute product T = ( Cgnu^T Cgmu )
    double complex T[ NSNS ] ;
    identity_spinmatrix( T ) ;
    gamma_spinmatrix( T , CgnuT ) ;
    spinmatrix_gamma( T , Cgmu ) ;

    // compute the trace of the above product
    double complex tr = spinmatrix_trace( T ) ;

    // loop open dirac indices
    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	// term[0] is proportional to the trace
	if( d1 == d2 ) {
	  mu_assert( "[UNIT] error : bar_contract_site term[0] broken" ,
		     cabs( term[0][ d2 + d1 * NS ] - ( 2.0 * NC * tr ) ) < FLTOL ) ;
	} else {
	  mu_assert( "[UNIT] error : bar_contract_site term[0] broken" ,
		     cabs( term[0][ d2 + d1 * NS ] ) < FLTOL ) ;
	}
	// term[1] is proportional to the product of the gamma matrices
	mu_assert( "[UNIT] error : bar_contract_site term[0] broken" ,
		   cabs( term[1][ d2 + d1 * NS ] - 2 * NC * T[ d2 + d1 * NS ] ) < FLTOL ) ;
	
      }
    }
  }

  // free the allocated memory
  free( term[1] ) ;
  free( term[0] ) ;
  free( term ) ;

  free( GAMMAS ) ;

  return NULL ;
}

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

  // check the higher level function
  mu_run_test( baryon_contract_site_test ) ;

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
