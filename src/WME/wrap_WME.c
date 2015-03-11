/**
   @file wrap_WME.c
   @brief wrapper for the WME calculations
 */

#include "common.h"

#include "WME.h"

// perform a WME contraction
int
contract_WME( FILE **fprops ,
	      const struct WME_info *wme ,
	      const int nWME )
{
  printf( "\n[WME] performing %d contraction(s) \n" , nWME ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nWME ; measurements++ ) {
    // loop the WME measurements
    if( WME( fprops[ wme[ measurements ].map[0] ] , 
	     wme[ measurements ].proptype1 ,
	     fprops[ wme[ measurements ].map[1] ] , 
	     wme[ measurements ].proptype2 ,
	     fprops[ wme[ measurements ].map[2] ] , 
	     wme[ measurements ].proptype3 ,
	     fprops[ wme[ measurements ].map[3] ] , 
	     wme[ measurements ].proptype4 ,
	     wme[ measurements ].outfile ) == FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
