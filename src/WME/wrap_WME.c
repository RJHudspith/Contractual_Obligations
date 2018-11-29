/**
   @file wrap_WME.c
   @brief wrapper for the WME calculations
 */
#include "common.h"

#include "GLU_timer.h"       // print_time()
#include "read_propheader.h" // reread the propagator file header
#include "WME.h"

// perform a WME contraction
int
contract_WME( struct propagator *prop ,
	      const struct WME_info *wme ,
	      const struct cut_info CUTINFO ,
	      const size_t nWME )
{
  fprintf( stdout , "\n[WME] performing %zu contraction(s) \n" , nWME ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nWME ; measurements++ ) {

    const size_t p1 = wme[ measurements ].map[0] ; // swall_0
    const size_t p2 = wme[ measurements ].map[1] ; // dwall_0
    const size_t p3 = wme[ measurements ].map[2] ; // swall_L/2
    const size_t p4 = wme[ measurements ].map[3] ; // dwall_L/2

    if( prop[p1].Source.type != WALL || prop[p2].Source.type != WALL ||
	prop[p3].Source.type != WALL || prop[p4].Source.type != WALL ) {
      fprintf( stderr , "[WME] Routine only works for WALL sources \n" ) ;
      return FAILURE ;
    }

    // loop the WME measurements
    if( WME( prop[p1] , prop[p2] , prop[p3] , prop[p4] ,
	     wme[ measurements ].outfile ) == FAILURE ) {
      return FAILURE ;
    }
    print_time() ;

    // reread headers
    if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
    if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
    if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
    if( reread_propheaders( &prop[ p4 ] ) == FAILURE ) { return FAILURE ; }
  }
  return SUCCESS ;
}
