/**
   @file wrap_WME.c
   @brief wrapper for the WME calculations
 */

#include "common.h"

#include "WME.h"

// perform a WME contraction
int
contract_WME( struct propagator *prop ,
	      const struct WME_info *wme ,
	      const int nWME )
{
  printf( "\n[WME] performing %d contraction(s) \n" , nWME ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nWME ; measurements++ ) {

    const int p0 = wme[ measurements ].map[0] ; // swall_0
    const int p1 = wme[ measurements ].map[1] ; // dwall_0
    const int p2 = wme[ measurements ].map[2] ; // swall_L/2
    const int p3 = wme[ measurements ].map[3] ; // dwall_L/2

    if( prop[p0].source != WALL || prop[p1].source != WALL ||
	prop[p2].source != WALL || prop[p3].source != WALL ) {
      printf( "[WME] routine only works for WALL sources \n" ) ;
      return FAILURE ;
    }

    // loop the WME measurements
    if( WME( prop[p0] , prop[p1] , prop[p2] , prop[p3] ,
	     wme[ measurements ].outfile ) == FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
