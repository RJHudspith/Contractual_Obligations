/**
   @file matops_tests.c
   @brief matrix operations tests
*/

#include "common.h"

#include "contractions.h"  // contractions
#include "gammas.h"
#include "minunit.h"       // minimal unit testing framework
#include "spinor_ops.h"    // identity_spinor()

#define FTOL ( NC * 1.E-14 ) 

static struct spinor A , Id ;         // single spinor storage

double complex *D = NULL ;

static struct gamma *GAMMAS = NULL ;

// test spincolor adjoint
static char 
*adjoint_test( void )
{
  struct spinor B ;
  adjoint_spinor( &B , A ) ;
  int d1 , d2 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  const double complex r = B.D[ d1 ][ d2 ].C[ c1 ][ c2 ] ;
	  const double complex s = A.D[ d2 ][ d1 ].C[ c2 ][ c1 ] ;
	  mu_assert("[CONTRACT UNIT] error : adjoint_spinor broken", 
		    !( fabs( creal( r ) - creal( s ) ) > FTOL ||
		       fabs( cimag( r ) + cimag( s ) ) > FTOL ) ) ;
	}
      }
    }
  }
  return NULL ;
}

// trace of a spinor Tr[ A * Id ] = Tr[ A ] 
static char *
bilinear_trace_test( void )
{
  const double complex tr = bilinear_trace( A , Id ) ;
  const double sol = NS * ( NC * ( NC - 1 ) * ( NC + 1 ) / 2 ) ;
  mu_assert("[CONTRACT UNIT] error : bilinear_trace broken", 
	    !( fabs( creal( tr ) - sol ) > FTOL ||
	       fabs( cimag( tr ) - sol ) > FTOL ) ) ;
  return NULL ;
}

// sanity check gamma matrices
static char 
*gammas_test( void )
{
  mu_assert("[CONTRACT UNIT] error : gamma matrix setup failure ", 
	    !( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) ) ;
  return NULL ;
}

// check full adjoint adj( gamma_5 Id gamma_5 ) = Id
static char *
full_adj_test( void )
{
  struct spinor adj ;
  full_adj( &adj , Id , GAMMAS[ GAMMA_5 ] ) ;
  const double *r = (const double*)adj.D ;
  const double *s = (const double*)Id.D ;
  size_t i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert("[CONTRACT UNIT] error : full_adj broken", 
	      !( fabs( *r - *s ) > FTOL ) ) ;
    r++ , s++ ;
  }
  return NULL ;
}

// test the left gamma multiply B = \gamma_testmat Id = \gamma_testmat
static char *
gamma_mul_l_test( void )
{
  const size_t testmat = GAMMA_1 ;
  struct spinor B = Id ;
  gamma_mul_l( &B , GAMMAS[ testmat ] ) ;
  size_t d1 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    const size_t col = (size_t)GAMMAS[ testmat ].ig[ d1 ] ;
    double complex comp = 1 ;
    switch( GAMMAS[ testmat ].g[ d1 ] ) {
    case 0 : comp =  1 ; break ;
    case 1 : comp =  I ; break ;
    case 2 : comp = -1 ; break ;
    case 3 : comp = -I ; break ;
    }
    for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	const double complex r = B.D[d1][col].C[c1][c2] ;
	if( c1 == c2 ) {
	  mu_assert("[CONTRACT UNIT] error : gamma_mul_l broken (diag)", 
		    !( fabs( creal( r ) - creal( comp ) ) > FTOL ||
		       fabs( cimag( r ) - cimag( comp ) ) > FTOL ) ) ;
	} else {
	  mu_assert("[CONTRACT UNIT] error : gamma_mul_l broken (off diag)", 
		    !( fabs( creal( r ) ) > FTOL ||
		       fabs( cimag( r ) ) > FTOL ) ) ;
	}
      }
    }
  }
  return NULL ;
}

// right multiply
static char *
gamma_mul_r_test( void )
{
  const int testmat = GAMMA_1 ;
  struct spinor B = Id ;
  gamma_mul_r( &B , GAMMAS[ testmat ] ) ;
  size_t d1 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    const size_t col = (size_t)GAMMAS[ testmat ].ig[ d1 ] ;
    double complex comp = 1 ;
    switch( GAMMAS[ testmat ].g[ col ] ) {
    case 0 : comp =  1 ; break ;
    case 1 : comp =  I ; break ;
    case 2 : comp = -1 ; break ;
    case 3 : comp = -I ; break ;
    }
    for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	const double complex r = B.D[col][d1].C[c1][c2] ;
	if( c1 == c2 ) {
	  mu_assert("[CONTRACT UNIT] error : gamma_mul_r broken (diag)", 
		    !( fabs( creal( r ) - creal( comp ) ) > FTOL ||
		       fabs( cimag( r ) - cimag( comp ) ) > FTOL ) ) ;
	} else {
	  mu_assert("[CONTRACT UNIT] error : gamma_mul_r broken (off diag)", 
		    !( fabs( creal( r ) ) > FTOL ||
		       fabs( cimag( r ) ) > FTOL ) ) ;
	}
      }
    }
  }
  return NULL ;
}

// test the left right multiply ( \gamma_l A \gamma_r ) 
// by using the independent left-right multiplies
static char *
gamma_mul_lr_test( void )
{
  struct spinor B = A ;
  gamma_mul_lr( &B , GAMMAS[ GAMMA_0 ] , GAMMAS[ GAMMA_1 ] ) ;

  struct spinor C = A ;
  gamma_mul_l( &C , GAMMAS[ GAMMA_0 ] ) ;
  gamma_mul_r( &C , GAMMAS[ GAMMA_1 ] ) ;

  const double *r = (const double*)B.D ;
  const double *s = (const double*)C.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert("[CONTRACT UNIT] error : gamma_mul_lr broken", 
	      !( fabs( *r - *s ) > FTOL ) ) ;
    r++ , s++ ;
  }
  return NULL ;
}

// compute the trace of a gamma matrix
static double complex
gamma_trace( const struct gamma G )
{
  size_t i ;
  double complex sum = 0.0 ;
  for( i = 0 ; i < NS ; i++ ) {
    if( G.ig[i] == i ) {
      switch( G.g[i] ) {
      case 0 : sum += +1 ; break ;
      case 1 : sum += +I ; break ;
      case 2 : sum += -1 ; break ;
      case 3 : sum += -I ; break ;
      }
    }
  }
  return sum ;
}

// test the meson contraction using our brute force methods
// contracts Tr[ \gamma_i ( \gamma_5 Id \gamma_5 )^\dagger \gamma_j Id
// with Tr[ \gamma_i \gamma_j ] for all gamma combinations
static char *
meson_contract_test( void )
{
  struct spinor adj ;
  const struct gamma G5 = GAMMAS[ GAMMA_5 ] ;
  size_t G1 , G2 ;
  for( G2 = 0 ; G2 < NSNS ; G2++ ) {
    for( G1 = 0 ; G1 < NSNS ; G1++ ) {

      // source and sink gammas
      struct gamma GSNK = GAMMAS[ G1 ] ;
      struct gamma GSRC = GAMMAS[ G2 ] ;

      // test against the product of the two inner gamma matrices
      struct gamma res ;
      gamma_mmul( &res , GSNK , GSRC ) ;
      double complex test = NC * gamma_trace( res ) ;

      // meson contract is
      // Tr[ GSNK ( G5 A G5 )^{\dagger} GSRC ( A ) ]
      const double complex tr1 = meson_contract( GSNK , Id , 
						 GSRC , Id , G5 ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !cabs( test - tr1 ) > FTOL ) ;

      full_adj( &adj , Id , G5 ) ;
      gamma_mul_l( &adj , GSNK ) ;
      gamma_mul_r( &adj , GSRC ) ;
      const double complex tr2 = bilinear_trace( adj , Id ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !( fabs( creal( tr1 ) - creal( tr2 ) ) > FTOL ||
		    fabs( cimag( tr1 ) - cimag( tr2 ) ) > FTOL ) ) ;

      full_adj( &adj , Id , G5 ) ;
      gamma_mul_lr( &adj , GSNK , GSRC ) ;
      const double complex tr3 = bilinear_trace( adj , Id ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !( fabs( creal( tr1 ) - creal( tr3 ) ) > FTOL ||
		    fabs( cimag( tr1 ) - cimag( tr3 ) ) > FTOL ) ) ;

      full_adj( &adj , Id , G5 ) ;
      const double complex tr4 = simple_meson_contract( GSNK , adj , 
							GSRC , Id ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !( fabs( creal( tr1 ) - creal( tr4 ) ) > FTOL ||
		    fabs( cimag( tr1 ) - cimag( tr4 ) ) > FTOL ) ) ;

    }
  }
  
  return NULL ;
}

// a more simple version of the above
static char *
simple_meson_contract_test( void )
{
  struct spinor adj ;
  struct gamma res ;
  size_t G1 , G2 ;
  for( G2 = 0 ; G2 < NSNS ; G2++ ) {
    for( G1 = 0 ; G1 < NSNS ; G1++ ) {

      // source and sink gammas
      const struct gamma GSNK = GAMMAS[ G1 ] ;
      const struct gamma GSRC = GAMMAS[ G2 ] ;

      gamma_mmul( &res , GSNK , GSRC ) ;
      double complex test = NC * gamma_trace( res ) ;

      // meson contract is
      // Tr[ GSNK A  GSRC ( A ) ]
      const double complex tr1 = simple_meson_contract( GSNK , Id , 
							GSRC , Id ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !cabs( test - tr1 ) > FTOL ) ;

      // left-right multiply
      adj = Id ;
      gamma_mul_l( &adj , GSNK ) ;
      gamma_mul_r( &adj , GSRC ) ;
      const double complex tr2 = bilinear_trace( adj , Id ) ;

      mu_assert( "[CONTRACT UNIT] error : meson_contract broken", 
		 !( fabs( creal( tr1 ) - creal( tr2 ) ) > FTOL ||
		    fabs( cimag( tr1 ) - cimag( tr2 ) ) > FTOL ) ) ;

      // left-right multiply
      adj = Id ;
      gamma_mul_lr( &adj , GSNK , GSRC ) ;
      const double complex tr3 = bilinear_trace( adj , Id ) ;

      mu_assert("[CONTRACT UNIT] error : simple_meson_contract broken", 
		!( fabs( creal( tr1 ) - creal( tr3 ) ) > FTOL ||
		   fabs( cimag( tr1 ) - cimag( tr3 ) ) > FTOL ) ) ;

    }
  }
  
  return NULL ;
}

// spinor tests
static char *
contractions_test( void )
{
  // initialise A as color matrices with values as their lexicographical
  // index
  double complex *b = (double complex*)A.D ;
  size_t i , j ;
  for( j = 0 ; j < NSNS ; j++ ) {
    for( i = 0 ; i < NCNC ; i++ ) {
      *b = ( i + I * i ) , b++ ;
    }
  }
  
  // set up an identity spinor
  identity_spinor( &Id ) ;

  // test gammas again, or for the first time depending on ordering
  mu_run_test( gammas_test ) ;
  mu_run_test( adjoint_test ) ;
  mu_run_test( bilinear_trace_test ) ;
  mu_run_test( gamma_mul_l_test ) ;
  mu_run_test( gamma_mul_r_test ) ;
  mu_run_test( gamma_mul_lr_test ) ;
  mu_run_test( full_adj_test ) ; // full adj uses gamma muls
                                 // test them first !!
  mu_run_test( simple_meson_contract_test ) ;
  mu_run_test( meson_contract_test ) ;

  return NULL ;
}

// driver for the tests
int
contractions_test_driver( void )
{
  // initialise the test counters
  tests_run = tests_fail = 0 ;

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // allocate spinmatrix "D"
  corr_malloc( (void**)&D , 16 , NSNS * sizeof( double complex ) ) ;

  // matrix operations test
  char *contres = contractions_test( ) ;

  free( D ) ;
  free( GAMMAS ) ;

  if( tests_fail == 0 ) {
    printf( "[CONTRACT UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , contres ) ;
    printf( "[CONTRACT UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
