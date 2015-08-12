/**
   @file bar_projections_tests.c
   @brief baryon_projection tests

   TODO :: full baryon projection unit test - J
 */
#include "common.h"

#include "bar_projections.h"  // projs, psqs
#include "gammas.h"           // make_gammas
#include "minunit.h"          // unit test framework
#include "spinmatrix_ops.h"   // spinmatrix operations

// gamma matrix technology
static struct gamma *GAMMA = NULL ; 

// temporary for pslash
static double complex *pslash = NULL ;
static double complex *t1 = NULL ;
static double complex *t2 = NULL ;
static double complex *t3 = NULL ;
static double complex *sum = NULL ;

// momenta stuff
static struct veclist *momentum = NULL ;
const size_t pidx = 0 ;

// write out a spinmatrix
static void
write_spinmatrix( double complex p[ NSNS ] )
{
  size_t d1 , d2 ;
  printf( "\n" ) ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      printf( " {%1.2f,%1.2f} " , 
	      creal( p[ d2 + d1*NS ] ) ,
	      cimag( p[ d2 + d1*NS ] ) ) ;
    }
    printf( "\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}

// momentum test
static char*
compute_p_psq_test( void )
{
  double p2 = 0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;
  mu_assert( "[UNIT] error : bar_projections compute p_psq\n" ,
	     fabs( p2 - (double)ND ) < 1E-14 ) ;
  return NULL ;
}

// driver for checking the various idempotency tests
static int
idempotent( void (*p1)( double complex *proj , const size_t i , const size_t j , 
			const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx ) , 
	    void (*p2)( double complex *proj , const size_t i , const size_t j , 
			const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx ) , 
	    void (*p3)( double complex *proj , const size_t i , const size_t j , 
			const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx ) , 
	    const size_t mu , 
	    const size_t nu ,
	    const struct gamma *GAMMA , 
	    const struct veclist *momentum , 
	    const size_t pidx )
{
  // set sum to zero
  zero_spinmatrix( sum ) ;
  // loop dummy index "rho"
  size_t rho ;
  for( rho = 0 ; rho < ND ; rho++ ) {
    p1( t1 , mu , rho , GAMMA , momentum , pidx ) ;
    p2( t2 , rho , nu , GAMMA , momentum , pidx ) ;
    spinmatrix_multiply( t3 , t1 , t2 ) ;
    atomic_add_spinmatrices( sum , t3 ) ;
  }
  p3( t1 , mu , nu , GAMMA , momentum , pidx ) ;
  // compare the product with what we expect it to be
  size_t d ;
  for( d = 0 ; d < NSNS ; d++ ) {
    if( cabs( t1[ d ] - sum[ d ] ) > 1E-14 ) {
      printf( "idempotency broken %zu %zu \n" , mu , nu ) ;
      write_spinmatrix( t3 ) ;
      write_spinmatrix( sum ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// ( P^{1/2}_{11} )_{ij} ( P^{J}_{cd} )_{jk} = 
//        \delta^{1/2,J} \delta^{1c} ( P^{1/2}_{ad} )_{ik}
static char*
p11_idempotency_test( void )
{
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      mu_assert( "[UNIT] error :: bar_projections P11.P11 != P11\n" , 
		 idempotent( P11 , P11 , P11 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P11.P12 != P12\n" , 
		 idempotent( P11 , P12 , P12 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P11.P21 != 0\n" , 
		 idempotent( P11 , P21 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P11.P22 != 0\n" , 
		 idempotent( P11 , P22 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P11.P32 != 0\n" , 
		 idempotent( P11 , P32 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
    }
  }
  return NULL ;
}

// ( P^{1/2}_{12} )_{ij} ( P^{J}_{cd} )_{jk} = 
//        \delta^{1/2,J} \delta^{2c} ( P^{1/2}_{ad} )_{ik}
static char*
p12_idempotency_test( void )
{
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      mu_assert( "[UNIT] error :: bar_projections P12.P11 != 0\n" , 
		 idempotent( P12 , P11 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P12.P12 != 0\n" , 
		 idempotent( P12 , P12 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P12.P21 != P11\n" , 
		 idempotent( P12 , P21 , P11 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P12.P22 != P12\n" , 
		 idempotent( P12 , P22 , P12 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P12.P32 != 0\n" , 
		 idempotent( P12 , P32 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
    }
  }
  return NULL ;
}

// ( P^{1/2}_{21} )_{ij} ( P^{J}_{cd} )_{jk} = 
//        \delta^{1/2,J} \delta^{1c} ( P^{1/2}_{ad} )_{ik}
static char*
p21_idempotency_test( void )
{
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      mu_assert( "[UNIT] error :: bar_projections P21.P11 != P21\n" , 
		 idempotent( P21 , P11 , P21 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P21.P12 != P22\n" , 
		 idempotent( P21 , P12 , P22 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P21.P21 != 0\n" , 
		 idempotent( P21 , P21 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P21.P22 != 0\n" , 
		 idempotent( P21 , P22 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P21.P32 != 0\n" , 
		 idempotent( P21 , P32 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
    }
  }
  return NULL ;
}

// ( P^{1/2}_{22} )_{ij} ( P^{J}_{cd} )_{jk} = 
//        \delta^{1/2,J} \delta^{2c} ( P^{1/2}_{ad} )_{ik}
static char*
p22_idempotency_test( void )
{
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      mu_assert( "[UNIT] error :: bar_projections P22.P11 != 0\n" , 
		 idempotent( P22 , P11 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P22.P12 != 0\n" , 
		 idempotent( P22 , P12 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P22.P21 != P21\n" , 
		 idempotent( P22 , P21 , P21 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P22.P22 != P22\n" , 
		 idempotent( P22 , P22 , P22 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P22.P32 != 0\n" , 
		 idempotent( P22 , P32 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
    }
  }
  return NULL ;
}

// ( P^{3/2}_{ab} )_{ij} ( P^{J}_{cd} )_{jk} = 
//        \delta^{3/2,J} \delta^{bc} ( P^{3/2}_{ad} )_{ik}
static char*
p32_idempotency_test( void )
{
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      mu_assert( "[UNIT] error :: bar_projections P32.P11 != 0\n" , 
		 idempotent( P32 , P11 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P32.P12 != 0\n" , 
		 idempotent( P32 , P12 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P32.P21 != 0\n" , 
		 idempotent( P32 , P21 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P32.P22 != 0\n" , 
		 idempotent( P32 , P22 , P00 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
      mu_assert( "[UNIT] error :: bar_projections P32.P32 != P32\n" , 
		 idempotent( P32 , P32 , P32 , mu , nu , 
			     GAMMA , momentum , pidx ) == SUCCESS ) ; 
    }
  }
  return NULL ;
}

static int
pslash_driver( void (*proj)( double complex *proj , const size_t i , 
			   const size_t j , const struct gamma *GAMMA , 
			   const struct veclist *momentum , const size_t pidx ) ,
	       double factor )
{
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;
  compute_pslash( pslash , GAMMA , p ) ;
  // loop all possible options
  size_t mu , nu , d ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      proj( t3 , mu , nu , GAMMA , momentum , 0 ) ;
      spinmatrix_multiply( t1 , pslash , t3 ) ;
      spinmatrix_multiply( t2 , t3 , pslash ) ;
      for( d = 0 ; d < NSNS ; d++ ) {
	if( cabs( t1[ d ] + factor * t2[ d ] ) > 1E-14 ) {
	  return FAILURE ;
	}
      }
    }
  }
  return SUCCESS ;
}

// check the pslash identities
static char *
check_pslashes( void )
{
  mu_assert( "[UNIT] error : bar_projections pslash.P11 != P11.pslash\n" ,
	     pslash_driver( P11 , -1 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections pslash.P12 != -P12.pslash\n" ,
	     pslash_driver( P12 , 1 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections pslash.P21 != -P12.pslash\n" ,
	     pslash_driver( P21 , 1 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections pslash.P22 != P22.pslash\n" ,
	     pslash_driver( P22 , -1 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections pslash.P32 != P32.pslash\n" ,
	     pslash_driver( P32 , -1 ) == SUCCESS ) ;
  return NULL ;
}

static int
pidentities_driver( void (*P)( double complex *proj , const size_t i , 
			       const size_t j , const struct gamma *GAMMA , 
			       const struct veclist *momentum , const size_t pidx ) )
{
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;
  size_t mu , nu , d ;
  for( nu = 0 ; nu < ND ; nu++ ) {
    zero_spinmatrix( sum ) ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      P( t1 , mu , nu , GAMMA , momentum , 0 ) ;
      spinmatrix_mulconst( t1 , p[ mu ] ) ;
      atomic_add_spinmatrices( sum , t1 ) ;
    }
    for( d = 0 ; d < ND ; d++ ) {
      if( cabs( sum[ d ] ) > 1E-14 ) {
	return FAILURE ;
      }
    }
  }
  return SUCCESS ;
}

// check p_\mu P = 0, for p11, p12 , p21 , p2
static char *
check_pidentities( void )
{
  mu_assert( "[UNIT] error : bar_projections p_mu P11 !=0\n" ,
	     pidentities_driver( P11 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections p_mu P12 !=0\n" ,
	     pidentities_driver( P12 ) == SUCCESS ) ;
  mu_assert( "[UNIT] error : bar_projections p_mu P11 !=0\n" ,
	     pidentities_driver( P32 ) == SUCCESS ) ;
  return NULL ;
}

// spinor tests
static char *
bar_projections_test( void )
{
  // check the momentum calc
  mu_run_test( compute_p_psq_test ) ;

  // run the idempotency tests
  mu_run_test( p11_idempotency_test ) ;
  mu_run_test( p12_idempotency_test ) ;
  mu_run_test( p21_idempotency_test ) ;
  mu_run_test( p22_idempotency_test ) ;
  mu_run_test( p32_idempotency_test ) ;

  // run the pslash tests
  mu_run_test( check_pslashes ) ;

  // run the momentum tests
  mu_run_test( check_pidentities ) ;

  return NULL ;
}

// runs the whole #!
int
bar_projections_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;
  
  // precompute the gamma basis
  GAMMA = malloc( NSNS * sizeof( struct gamma ) ) ;
  corr_malloc( (void**)&pslash , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&t1, 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&t2 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&t3 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&sum , 16 , NSNS * sizeof( double complex ) ) ;

  // make the gammas in the CHIRAL basis
  make_gammas( GAMMA , CHIRAL ) ;

  // set up veclist as momenta of all ones
  momentum = malloc( sizeof( struct veclist ) ) ;
  momentum[0].idx = 0 ;
  momentum[0].nsq = 0 ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    momentum[ 0 ].MOM[ mu ] = 1 ;
    momentum[ 0 ].nsq += 1 ;
  }

  // initial gamma setup and test
  char *barprojres = bar_projections_test( ) ;

  // memfree
  free( momentum ) ;
  free( sum ) ;
  free( t3 ) ;
  free( t2 ) ;
  free( t1 ) ;
  free( pslash ) ;
  free( GAMMA ) ;

  if( tests_fail == 0 ) {
    printf( "[BAR PROJECTIONS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , barprojres ) ;
    printf( "[BAR PROJECTIONS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
