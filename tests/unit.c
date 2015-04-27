/**
   @file unit.c
   @brief very simple unit testing framework
 */

#include "common.h"

#include "geometry.h"      // init_geom()
#include "spinor_tests.h"

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
  // have a look at spinor operations
  spinor_test_driver( ) ;

  return SUCCESS ;
}
