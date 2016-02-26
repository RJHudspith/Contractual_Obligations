/**
   @file wrap_WME.c
   @brief wrapper for the WME calculations
 */
#include "common.h"

#include "WME.h"
#include "read_propheader.h" // reread the propagator file header

// perform a WME contraction
int
contract_WME( struct propagator *prop ,
	      const struct WME_info *wme ,
	      const struct cut_info CUTINFO ,
	      const size_t nWME )
{
  printf( "\n[WME] performing %zu contraction(s) \n" , nWME ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nWME ; measurements++ ) {

    const size_t p0 = wme[ measurements ].map[0] ; // swall_0
    const size_t p1 = wme[ measurements ].map[1] ; // dwall_0
    const size_t p2 = wme[ measurements ].map[2] ; // swall_L/2
    const size_t p3 = wme[ measurements ].map[3] ; // dwall_L/2

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
    // reread headers 
    rewind( prop[p0].file ) ; read_propheader( &prop[p0] ) ;
    rewind( prop[p1].file ) ; read_propheader( &prop[p1] ) ;
    rewind( prop[p2].file ) ; read_propheader( &prop[p2] ) ;
    rewind( prop[p3].file ) ; read_propheader( &prop[p3] ) ;
  }
  return SUCCESS ;
}
