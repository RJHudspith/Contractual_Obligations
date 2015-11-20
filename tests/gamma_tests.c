#include "common.h"

#include "gammas.h"
#include "minunit.h"

static struct gamma *GAMMAS = NULL ; // gamma matrix technology

// gamma matrix comparison
static int
gamma_comparison( const struct gamma G1 , 
		  const struct gamma G2 )
{
  int j ;
  for( j = 0 ; j < NS ; j++ ) {
    if( G1.ig[j] != G2.ig[j] || G1.g[j] != G2.g[j] ) {
      return 1 ;
    }
  }
  return 0 ;
}

// check the multiplies
static char *
gamma_mmul_test( void )
{
  struct gamma res ;
  int mu ;
  for( mu = 0 ; mu < NSNS ; mu++ ) {
    gamma_mmul( &res , GAMMAS[ IDENTITY ] , GAMMAS[ mu ] ) ;
    mu_assert( "[UNIT] error : gammas gamma_mmul broken on identity multiply1" ,
	       !gamma_comparison( res , GAMMAS[ mu ] ) ) ;
    gamma_mmul( &res , GAMMAS[ mu ] , GAMMAS[ IDENTITY ] ) ;
    mu_assert( "[UNIT] error : gammas gamma_mmul broken on identity multiply2" ,
	       !gamma_comparison( res , GAMMAS[ mu ] ) ) ;
  }
  // check the gammas square to the identity
  gamma_mmul( &res , GAMMAS[ GAMMA_0 ] , GAMMAS[ GAMMA_0 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g0.g0 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_1 ] , GAMMAS[ GAMMA_1 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g1.g1 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_2 ] , GAMMAS[ GAMMA_2 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g2.g2 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_3 ] , GAMMAS[ GAMMA_3 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g3.g3 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_5 ] , GAMMAS[ GAMMA_5 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g5.g5 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;

  // compute gamma_0.gamma_1.gamma_2.gamma_3 = gamma_5 and compare
  struct gamma t1 , t2 ;
  gamma_mmul( &t1 , GAMMAS[ GAMMA_0 ] , GAMMAS[ GAMMA_1 ] ) ;
  gamma_mmul( &t2 , GAMMAS[ GAMMA_2 ] , GAMMAS[ GAMMA_3 ] ) ;
  gamma_mmul( &res , t1 , t2 ) ;
  mu_assert( "[UNIT] error : g0.g1.g2.g3 != g5" ,
	     !gamma_comparison( res , GAMMAS[ GAMMA_5 ] ) ) ;
  return NULL ;
}

// test if we can take a conjugate
static char*
gconj_test( void )
{
  mu_assert( "[UNIT] error : gconj 0" , ( gconj( 0 ) == 0 ) ) ;
  mu_assert( "[UNIT] error : gconj 1" , ( gconj( 1 ) == 3 ) ) ;
  mu_assert( "[UNIT] error : gconj 2" , ( gconj( 2 ) == 2 ) ) ;
  mu_assert( "[UNIT] error : gconj 3" , ( gconj( 3 ) == 1 ) ) ;
  return NULL ;
}

// test that conjugate of identity is the identity
static char*
gamma_conj_test( void )
{
  mu_assert( "[UNIT] error : gammas (I*).I != I" ,
	     !( gamma_comparison( gamma_conj( GAMMAS[ IDENTITY ] ) ,
				  GAMMAS[ IDENTITY ] ) ) ) ;
  return NULL ;
}

// test that the dagger of a gamma matrix is the gamma matrix back
static char*
gamma_dag_test( void )
{
  size_t mu ;
  for( mu = 0 ; mu < GAMMA_5 ; mu++ ) {
    mu_assert( "[UNIT] error : g^{dagger} != g" , 
	       !( gamma_comparison( gamma_dag( GAMMAS[mu] ) ,
				    GAMMAS[mu] ) ) ) ;
  }
  return NULL ;
}

// test that the transpose of the identity is the identity
static char*
gamma_transpose_test( void )
{
  struct gamma IT = gamma_transpose( GAMMAS[ IDENTITY ] ) ;
  mu_assert( "[UNIT] error : I^T != I" , 
	     !( gamma_comparison( IT , GAMMAS[IDENTITY] ) ) ) ;
  return NULL ;
}

// test i^4*I
static char*
gamma_muli_test( void )
{
  struct gamma iI = GAMMAS[ IDENTITY ] ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  mu_assert( "[UNIT] error : gammas multiply by i" ,
	     !( gamma_comparison( iI , GAMMAS[ IDENTITY ] ) ) ) ;
  return NULL ;
}

// test Cgamma_mu = g1.g3 , according to gattringer and Lang C has
// the following properties C = C^-1 = C^{\dagger}
static char*
Cgmu_test( void )
{
  struct gamma C = CGmu( GAMMAS[ IDENTITY ] , GAMMAS ) ;
  // C = C^{dagger} ( A.25 )
  struct gamma Cdag = gamma_dag( C ) ;
  mu_assert( "[UNIT] error : C matrix is not own dagger" ,
	     !( gamma_comparison( C , Cdag ) ) ) ;
  // C = -C^{T} ( A.25 )
  struct gamma Ctrans = gamma_transpose( C ) ;
  gamma_muli( &Ctrans ) ;
  gamma_muli( &Ctrans ) ;
  mu_assert( "[UNIT] error : C matrix is not own dagger" ,
	     !( gamma_comparison( C , Ctrans ) ) ) ;
  // C \gamma_\mu C = -\gamma_\mu^T ( A.23 )
  size_t mu ;
  for( mu = 0 ; mu < IDENTITY ; mu++ ) {
    struct gamma t1 , t2 ;
    gamma_mmul( &t1 , C , GAMMAS[mu] ) ;
    gamma_mmul( &t2 , t1 , C ) ;
    struct gamma t3 = gamma_transpose( GAMMAS[mu] ) ;
    gamma_muli( &t3 ) ;
    gamma_muli( &t3 ) ;
    mu_assert( "[UNIT] error : Cg_mu C != g_mu^T" ,
	       !( gamma_comparison( t2 , t3 ) ) ) ;
  }
  return NULL ;
}

// test \gamma_t.G*.\gamma_t
static char*
gt_Gconj_gt_test( void )
{
  // test that with G==I we get back I
  mu_assert( "[UNIT] error : gt.I.gt != I" ,
	     !( gamma_comparison( GAMMAS[IDENTITY] ,
				  gt_Gconj_gt( GAMMAS[IDENTITY] ,
					       GAMMAS ) ) ) ) ;
  /* this is untrue!!!
  struct gamma Cg5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma gtCg5gt = gt_Gconj_gt( Cg5 , GAMMAS ) ;
  mu_assert( "[UNIT] error : gt.Cg5.gt != Cg5" ,
	     !( gamma_comparison( Cg5 , gtCg5gt ) ) ) ;
  */
  return NULL ;
}

// test \gamma_t.G*.\gamma_t
static char*
gt_Gdag_gt_test( void )
{
  // test that with G==I we get back I
  size_t mu ;
  for( mu = 0 ; mu < 6 ; mu++ ) {
    struct gamma res = GAMMAS[mu] ;
    // these two commute with gamma_3
    if( mu != GAMMA_3 && mu != IDENTITY) {
      gamma_muli( &res ) ;
      gamma_muli( &res ) ;
    }
    mu_assert( "[UNIT] error : gt.g_mu.gt != -g_mu" ,
	       !( gamma_comparison( res ,
				    gt_Gdag_gt( GAMMAS[mu] ,
						GAMMAS ) ) ) ) ;
  }
  //
  struct gamma Cg5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma gtCg5gt = gt_Gdag_gt( Cg5 , GAMMAS ) ;
  mu_assert( "[UNIT] error : gt.Cg5.gt != Cg5" ,
	     !( gamma_comparison( Cg5 , gtCg5gt ) ) ) ;
  return NULL ;
}

static char*
gammas_test( void )
{
  // independent tests
  mu_run_test( gamma_mmul_test ) ;
  mu_run_test( gconj_test ) ;
  mu_run_test( gamma_conj_test ) ;
  mu_run_test( gamma_dag_test ) ;
  mu_run_test( gamma_muli_test ) ;
  mu_run_test( gamma_transpose_test ) ;
  
  // these rely on the above
  mu_run_test( Cgmu_test ) ;
  mu_run_test( gt_Gconj_gt_test ) ;
  mu_run_test( gt_Gdag_gt_test ) ;

  return NULL ;
}

int
gamma_test_driver( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;
  
  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  make_gammas( GAMMAS , CHIRAL ) ;

  char *gammares = gammas_test( ) ;

  make_gammas( GAMMAS , NREL ) ;

  gammares = gammas_test( ) ;

  free( GAMMAS ) ;

  if( tests_fail == 0 ) {
    printf( "[GAMMAS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return SUCCESS ;
  } else {
    printf( "%s \n" , gammares ) ;
    printf( "[GAMMAS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }
}
