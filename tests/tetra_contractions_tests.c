/**
   @file tetra_contractions_tests.c
   @brief tetraquark contraction tests
 */
#include "common.h"

#include "contractions.h"       // simple_meson_contract
#include "gammas.h"             // Cgmu, make_gammas
#include "minunit.h"            // mu_assert
#include "spinor_ops.h"         // spinor_identity
#include "spinmatrix_ops.h"     // get_spinmatrix
#include "tetra_contractions.h" // precompute_block, get_abcd ...

// our tolerance
#define FLTOL (NC*1.E-14)

// some temporary space
static struct spinor S1 , S2 ;

// gamma space
static struct gamma *GAMMAS = NULL ;

// block matrices space
struct block *C1 = NULL ;

// compare a spinmatrix to a gamma matrix
static int
compare_spinmatrix_to_gamma( double complex *s ,
			     const struct gamma G )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      if( G.ig[ d1 ] == d2 ) {
	switch( G.g[d1] ) {
	case 0 : if( cabs( s[d2+NS*d1] - 1 ) > FLTOL ) ; return 1 ;
	case 1 : if( cabs( s[d2+NS*d1] - I ) > FLTOL ) ; return 1 ;
	case 2 : if( cabs( s[d2+NS*d1] + 1 ) > FLTOL ) ; return 1 ;
	case 3 : if( cabs( s[d2+NS*d1] + I ) > FLTOL ) ; return 1 ;
	}
      } else {
	if( cabs( s[d2+NS*d1] ) > FLTOL ) ; return 1 ;
      }
    }
  }
  return 0 ;
}

// compute the trace of a gamma object
static double complex
gamma_trace( const struct gamma G )
{
  register double complex sum = 0.0 ;
  size_t d1 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    if( G.ig[ d1 ] == d1 ) {
      switch( G.g[d1] ) {
      case 0 : sum += 1 ; break ;
      case 1 : sum += I ; break ;
      case 2 : sum -= 1 ; break ;
      case 3 : sum -= I ; break ;
      }
    }
  }
  return sum ;
}

// test our linearised index works
static char*
element_test( void )
{
  mu_assert( "[UNIT] error : element broken\n" ,
	     element( NC-1 , NC-1 , NC-1 , NC-1 ) == (NCNC*NCNC-1) ) ;
  return NULL ;
}

// test that we get the correct color elements
static char*
get_abcd_test( void )
{
  // test it is an invertible map
  size_t a , b , c , d , abcd ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    mu_assert( "[UNIT] error : get_abcd broken\n" ,
	       element( a , b , c , d ) == abcd ) ; 
  }
  return NULL ;
}

// put two identity spinors in the precompute block
static char*
precompute_block_test( void )
{
  size_t GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < B_CHANNELS ; GSNK++ ) {

      identity_spinor( &S1 ) ;
      identity_spinor( &S2 ) ;
    
      // common use case G1 == C_gamma_i
      // G2 == tilde{C}_gamma_5
      struct gamma Cgi = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
      struct gamma Cgj = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
      struct gamma tildeCgj = gt_Gdag_gt( Cgj , GAMMAS[ GAMMA_T] ) ;

      precompute_block( C1 , S1 , Cgi , S2 , tildeCgj ) ;
  
      // result should be spinmatrices that are the gamma product on 
      // the diagonal where a=b=c=d
      struct gamma prod ;
      gamma_mmul( &prod , Cgi , tildeCgj ) ;
  
      size_t abcd , a , b , c , d , d1d2 ;
      for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
	get_abcd( &a , &b , &c , &d , abcd ) ;
	// if we hit the identity components we compare to
	// the product of gammas
	if( a == b && c == d ) {
	  mu_assert( "[UNIT] error : precompute_block broken\n" ,
		     compare_spinmatrix_to_gamma( C1[abcd].M , prod ) ) ;
	} else {
	  // should be zero
	  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
	    mu_assert( "[UNIT] error : precompute_block broken\n" ,
		       cabs( C1[abcd].M[d1d2] ) < FLTOL ) ;
	  }
	}
	//
      }
    }
  }
  return NULL ;
}

// test that sum_{c,d} Tr_S[ Block[dccd].M ] is the same as a spin-color trace
static char *
spincolor_trace_test( void )
{
  size_t GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < B_CHANNELS ; GSNK++ ) {

      identity_spinor( &S1 ) ;
      identity_spinor( &S2 ) ;
    
      // common use case G1 == C_gamma_i
      // G2 == tilde{C}_gamma_5
      struct gamma Cgi = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
      struct gamma Cgj = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
      struct gamma tildeCgj = gt_Gdag_gt( Cgj , GAMMAS[ GAMMA_T ] ) ;
      precompute_block( C1 , S1 , Cgi , S2 , tildeCgj ) ;

      // product of gamma matrices
      struct gamma prod ;
      gamma_mmul( &prod , Cgi , tildeCgj ) ;
      double complex trprod = NC * gamma_trace( prod ) ;

      // perform explicit spin trace
      register double complex sum1 = 0.0 ;
      size_t c , d ;
      for( c = 0 ; c < NC ; c++ ) {
	for( d = 0 ; d < NC ; d++ ) {
	  sum1 += spinmatrix_trace( C1[ element( d , c , c , d ) ].M ) ;
	}
      }
      double complex sum2 = -simple_meson_contract( tildeCgj , S1 , Cgi , S2 ) ;

      mu_assert( "[UNIT] error : spincolor_trace_test broken #1\n" ,
		 cabs( sum1 - trprod ) < FLTOL ) ;
      mu_assert( "[UNIT] error : spincolor_trace_test broken #2\n" ,
		 cabs( sum1 - sum2 ) < FLTOL ) ;

      // for non-SU(3) the index structure is different and so this
      // test fails - J
#if NC == 3
      // test the spintrace-spintrace identity this is the first term in 
      // the mesonic contraction
      size_t abcd , a , b ;
      sum1 = sum2 = 0.0 ;
      for( abcd = 0 ; abcd < NSNS*NSNS ; abcd++ ) {
	get_abcd( &a , &b , &c , &d , abcd ) ;
	sum1 += 
	  spinmatrix_trace( C1[ element( d , a , a , d ) ].M ) * 
	  spinmatrix_trace( C1[ element( c , b , b , c ) ].M ) ;
      }
      // factors of 2 are because we double count with the 
      // full sum over a,b,c,d
      sum2 = 
	2 * simple_meson_contract( tildeCgj , S1 , Cgi , S2 ) * 
	simple_meson_contract( tildeCgj , S1 , Cgi , S2 ) ;
      trprod *= 2 * trprod ;  

      mu_assert( "[UNIT] error : spincolor_trace_test broken \n" ,
		 cabs( sum1 - trprod ) < FLTOL ) ;
      mu_assert( "[UNIT] error : spincolor_trace_test broken \n" ,
		 cabs( sum1 - sum2 ) < FLTOL ) ;
#endif
      //
    }
  }
  return NULL ;
}

// test that 
// S^{\alpha\beta}_{aa} I^{\beta\kappa} S^{\kappa\delta}_{aa} I^{\delta\alpha}
// is equal to the spin trace Tr_S[ S_{aa} G S_{aa} G ]
static char *
spinmatrix_trace_test( void )
{
  identity_spinor( &S1 ) ;
  identity_spinor( &S2 ) ;

  double complex s1[ NSNS ] , s2[ NSNS ] ;
  get_spinmatrix( s1 , S1 , 0 , 0 ) ;
  get_spinmatrix( s2 , S2 , 0 , 0 ) ;

  precompute_block( C1 , 
		    S1 , GAMMAS[ IDENTITY ] , 
		    S2 , GAMMAS[ IDENTITY ] ) ;
  size_t alpha , beta ;
  double complex sum1 = 0.0 , sum2 = 0.0 ;
  for( alpha = 0 ; alpha < NS ; alpha++ ) {
    for( beta = 0 ; beta < NS ; beta++ ) {
      sum1 += 
	s1[ beta + alpha * NS ] * 
	s2[ beta + alpha * NS ] ;
    }
  }
  sum2 = spinmatrix_trace( C1[ 0 ].M ) ;
  
  mu_assert( "[UNIT] error : spinmatrix_trace identity broken\n" ,
	     cabs( sum1 - sum2 ) < FLTOL ) ;

  return NULL ;
}

// baryon operations tests
static char *
tetra_contractions_test( void )
{
  // test our linearised tetraquark indices
  mu_run_test( element_test ) ;
  mu_run_test( get_abcd_test ) ;

  // check the higher level function
  mu_run_test( precompute_block_test ) ;
  mu_run_test( spincolor_trace_test ) ;
  mu_run_test( spinmatrix_trace_test ) ;

  return NULL ;
}

// runs the whole #!
int
tetra_contractions_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  make_gammas( GAMMAS , CHIRAL ) ;

  // allocate block
  corr_malloc( (void**)&C1 , 16 , NCNC*NCNC*sizeof( struct block ) ) ;

  // initial gamma setup and test
  char *tetrasres = tetra_contractions_test( ) ;

  free( C1 ) ;
  free( GAMMAS ) ;

  if( tests_fail == 0 ) {
    printf( "[TETRAS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , tetrasres ) ;
    printf( "[TETRAS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
