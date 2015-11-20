/**
   @file unit.c
   @brief very simple unit testing framework
 */

#include "common.h"

#include "bar_projections_tests.h"
#include "bar_ops_tests.h"
#include "contract_tests.h"  
#include "geometry.h"        // init_geom()
#include "matops_tests.h"
#include "spinmatrix_tests.h"
#include "spinor_tests.h"
#include "gamma_tests.h"

struct latt_info Latt ;

int tests_fail = 0 ;
int tests_run = 0 ;

int 
main( const int argc , const char *argv[] )
{
  // set up a very small geometry 4^{ND-1} (LCU)
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 4 ;
  }
  init_geom( ) ; // geometry inherited from GLU and has been tested

  printf( "\n" ) ;

  int total = 0 ;

  // have a look at the gamma operations
  if( gamma_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at spinor operations
  if( contractions_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at matrix operations
  if( matops_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at spinmatrix operations
  if( spinmatrix_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at spinor operations
  if( spinor_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at baryon projection operations
  if( bar_projections_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  // have a look at baryon operations
  if( bar_ops_test_driver( ) == FAILURE ) goto failure ;
  total += tests_run ;

  printf( "[UNIT] %d tests run and passed\n" , total ) ;

  return SUCCESS ;

 failure :
  return FAILURE ;
}
