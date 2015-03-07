/**
   @file wrap_mesons.c
   @brief meson contractions wrappers
 */

#include "common.h"

#include "brutal_mesons.h"   // brutal mesons
#include "mesons.h"          // meson contractions

// meson contraction driver
int
contract_mesons( FILE **fprops ,
		 const struct meson_info *mesons ,
		 const int nmesons )
{
  printf( "[MESONS] performing %d contractions \n" , nmesons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nmesons ; measurements++ ) {
    switch( mesons[ measurements ].source ) {
    case WALL :
      if( mesons[ measurements ].map[0] == mesons[ measurements ].map[0] ) {
	return wall_mesons( fprops[ mesons[ measurements ].map[0] ] , 
			    mesons[ measurements ].proptype1 ,
			    mesons[ measurements ].outfile ) ;
      } else {
	return wall_double_mesons( fprops[ mesons[ measurements ].map[0] ] , 
				   mesons[ measurements ].proptype1 ,
				   fprops[ mesons[ measurements ].map[1] ] , 
				   mesons[ measurements ].proptype2 ,
				   mesons[ measurements ].outfile ) ;
      }
      break ;
    case POINT :
      if( mesons[ measurements ].map[0] == mesons[ measurements ].map[0] ) {
	return single_mesons( fprops[ mesons[ measurements ].map[0] ] , 
			      mesons[ measurements ].proptype1 ,
			      mesons[ measurements ].outfile ) ;
      } else {
	return double_mesons( fprops[ mesons[ measurements ].map[0] ] , 
			      mesons[ measurements ].proptype1 ,
			      fprops[ mesons[ measurements ].map[1] ] , 
			      mesons[ measurements ].proptype2 ,
			      mesons[ measurements ].outfile ) ;
      }
      break ;
    }
  }
  return SUCCESS ;
}
