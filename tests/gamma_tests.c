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
  gamma_mmul( &res , GAMMAS[ GAMMA_X ] , GAMMAS[ GAMMA_X ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g0.g0 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_Y ] , GAMMAS[ GAMMA_Y ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g1.g1 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_Z ] , GAMMAS[ GAMMA_Z ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g2.g2 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_T ] , GAMMAS[ GAMMA_T ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g3.g3 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_5 ] , GAMMAS[ GAMMA_5 ] ) ;
  mu_assert( "[UNIT] error : gammas gamma_mmul broken g5.g5 != 1" ,
	     !gamma_comparison( res , GAMMAS[ IDENTITY ] ) ) ;

  // compute gamma_0.gamma_1.gamma_2.gamma_3 = gamma_5 and compare
  struct gamma t1 , t2 ;
  gamma_mmul( &t1 , GAMMAS[ GAMMA_X ] , GAMMAS[ GAMMA_Y ] ) ;
  gamma_mmul( &t2 , GAMMAS[ GAMMA_Z ] , GAMMAS[ GAMMA_T ] ) ;
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
  // test that gamma_conj is gamma dag transposed for all gammas
  size_t mu ;
  for( mu = 0 ; mu < GAMMA_5 ; mu++ ) {
    mu_assert( "[UNIT] error : gamma_conj (g*) != (g^{dagger})^T" ,
	     !( gamma_comparison( gamma_conj( GAMMAS[ mu ] ) ,
				  gamma_transpose( gamma_dag( GAMMAS[ mu ] ) ) 
				  ) ) ) ;
    
  }
  return NULL ;
}

// test that the dagger of a gamma matrix is the same gamma matrix
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

// transpose test
static char*
gamma_transpose_test( void )
{
  // test that the transpose of the identity is the identity
  struct gamma IT = gamma_transpose( GAMMAS[ IDENTITY ] ) ;
  mu_assert( "[UNIT] error : I^T != I" , 
	     !( gamma_comparison( IT , GAMMAS[IDENTITY] ) ) ) ;
  // test that ( gamma_i^T gamma_j^T ) = ( gamma_j gamma_i )^T
  size_t i , j ;
  struct gamma res1 , res2 ;
  for( i = 0 ; i < GAMMA_5 ; i++ ) {
    for( j = 0 ; j < GAMMA_5 ; j++ ) {
      gamma_mmul( &res1 , gamma_transpose( GAMMAS[ i ] ) ,
		  gamma_transpose( GAMMAS[ j ] ) ) ;
      gamma_mmul( &res2 , GAMMAS[ j ] , GAMMAS[ i ] ) ;
      mu_assert( "[UNIT] error : g_i^T.g_j^T != (g_2 g_1)^T" , 
		 !( gamma_comparison( res1 , gamma_transpose( res2 ) ) ) ) ;
    }
  }
  return NULL ;
}

// test i^4*I = I
static char*
gamma_muli_test( void )
{
  struct gamma iI = GAMMAS[ IDENTITY ] ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  gamma_muli( &iI ) ;
  mu_assert( "[UNIT] error : gamma_muli (i)^4 * I != I " ,
	     !( gamma_comparison( iI , GAMMAS[ IDENTITY ] ) ) ) ;
  return NULL ;
}

// test i^3*I = -i*I
static char*
gamma_mul_minus1_test( void )
{
  struct gamma iI1 = GAMMAS[ IDENTITY ] , iI2 = GAMMAS[ IDENTITY ] ;
  gamma_mul_minusi( &iI1 ) ;
  // i^3 = -i
  gamma_muli( &iI2 ) ;
  gamma_muli( &iI2 ) ;
  gamma_muli( &iI2 ) ;
  mu_assert( "[UNIT] error : gamma_muli (i)^3 * I != -i * I" ,
	     !( gamma_comparison( iI1 , iI2 ) ) ) ;
  return NULL ;
}

// test i^2*I = -I
static char*
gamma_mul_minusi_test( void )
{
  struct gamma iI1 = GAMMAS[ IDENTITY ] , iI2 = GAMMAS[ IDENTITY ] ;
  gamma_mul_minus1( &iI1 ) ;
  // i^3 = -i
  gamma_muli( &iI2 ) ;
  gamma_muli( &iI2 ) ;
  mu_assert( "[UNIT] error : gamma_muli (i)^2 * I != -I" ,
	     !( gamma_comparison( iI1 , iI2 ) ) ) ;
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
  mu_assert( "[UNIT] error : C != C^dagger" ,
	     !( gamma_comparison( C , Cdag ) ) ) ;
  // C = -C^{T} ( A.25 )
  struct gamma Ctrans = gamma_transpose( C ) ;
  gamma_mul_minus1( &Ctrans ) ;
  mu_assert( "[UNIT] error : C != -C^T" ,
	     !( gamma_comparison( C , Ctrans ) ) ) ;
  // C \gamma_\mu C = -\gamma_\mu^T ( A.23 )
  size_t mu ;
  for( mu = 0 ; mu < IDENTITY ; mu++ ) {
    struct gamma t1 , t2 ;
    gamma_mmul( &t1 , C , GAMMAS[mu] ) ;
    gamma_mmul( &t2 , t1 , C ) ;
    struct gamma t3 = gamma_transpose( GAMMAS[mu] ) ;
    gamma_mul_minus1( &t3 ) ;
    mu_assert( "[UNIT] error : Cg_mu C != g_mu^T" ,
	       !( gamma_comparison( t2 , t3 ) ) ) ;
  }
  // C.C = Id
  struct gamma CC ;
  gamma_mmul( &CC , C , C ) ;
  mu_assert( "[UNIT] error : C.C != I" ,
	     !( gamma_comparison( CC , GAMMAS[ IDENTITY ] ) ) ) ;
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
					       GAMMAS[GAMMA_T] ) ) ) ) ;
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
    // these two commute with gamma_t
    if( mu != GAMMA_T && mu != IDENTITY) {
      gamma_muli( &res ) ;
      gamma_muli( &res ) ;
    }
    mu_assert( "[UNIT] error : gt.g_mu.gt != -g_mu" ,
	       !( gamma_comparison( res ,
				    gt_Gdag_gt( GAMMAS[ mu ] ,
						GAMMAS[ GAMMA_T ] ) ) ) ) ;
  }
  // test this identity
  for( mu = 0 ; mu < 6 ; mu++ ) {
    if( mu == GAMMA_Y || mu == IDENTITY ) continue ;
    struct gamma Cgi = CGmu( GAMMAS[ mu ] , GAMMAS ) ;
    struct gamma gtCgigt = gt_Gdag_gt( Cgi , GAMMAS[ GAMMA_T ] ) ;
    mu_assert( "[UNIT] error : gt.Cgi.gt != Cgi" ,
	       !( gamma_comparison( Cgi , gtCgigt ) ) ) ;
  }
  // test that gt.(C.gy)^dagger.gt = igt
  struct gamma Cgy = CGmu( GAMMAS[ GAMMA_Y ] , GAMMAS ) ;
  struct gamma gtCgygt = gt_Gdag_gt( Cgy , GAMMAS[ GAMMA_T ] ) ;
  struct gamma igt = GAMMAS[ GAMMA_T ] ;
  gamma_muli( &igt ) ;
  mu_assert( "[UNIT] error : gt.Cgy.gt != igt" ,
	     !( gamma_comparison( gtCgygt , igt ) ) ) ;
  // test that gt.(C)^dagger.gt = -C
  struct gamma CgI = CGmu( GAMMAS[ IDENTITY ] , GAMMAS ) ;
  struct gamma gtCgIgt = gt_Gdag_gt( CgI , GAMMAS[ GAMMA_T ] ) ;
  mu_assert( "[UNIT] error : gt.CgI.gt != -C" ,
	     !( gamma_comparison( gtCgIgt , gamma_transpose( CgI ) ) ) ) ;
  // test that gt (C.gamma)^\dagger gt = gt gamma.C gt ?
  for( mu = 0 ; mu < 6 ; mu++ ) {
    struct gamma gt_Cgmudag_gt = gt_Gdag_gt( CGmu( GAMMAS[ mu ] , GAMMAS ) , 
					     GAMMAS[ GAMMA_T ] ) ;
    struct gamma tmp , res ;
    gamma_mmul( &res , CgI , GAMMAS[ GAMMA_T ] ) ;
    gamma_mmul( &tmp , GAMMAS[ mu ] , res ) ;
    gamma_mmul( &res , GAMMAS[ GAMMA_T ] , tmp ) ;
    mu_assert( "[UNIT] error : gt.(C.gmu)^dagger.gt != gt.(gmu.C).gt \n" ,
	       !( gamma_comparison( res , gt_Cgmudag_gt ) ) ) ;
  }
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
  mu_run_test( gamma_mul_minus1_test ) ;
  mu_run_test( gamma_mul_minusi_test ) ;
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

  if( tests_fail != 0 ) {
    printf( "%s \n" , gammares ) ;
    printf( "[GAMMAS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return FAILURE ;
  }

  make_gammas( GAMMAS , NREL_FWD ) ;

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
