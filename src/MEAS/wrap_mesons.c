/**
   @file wrap_mesons.c
   @brief meson contractions wrappers
 */

#include "common.h"

#include "brutal_mesons.h"   // brutal mesons
#include "wall_mesons.h"     // wall mesons is now the bona-fide mesons

// meson contraction driver
int
contract_mesons( struct propagator *prop ,
		 const struct meson_info *mesons ,
		 const int nmesons )
{
  printf( "\n[MESONS] performing %d contraction(s) \n" , nmesons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nmesons ; measurements++ ) {
    // to make it more legible
    const int p1 = mesons[ measurements ].map[0] ;
    const int p2 = mesons[ measurements ].map[1] ;

    if( p1 == p2 ) {
      if( mesons_diagonal( prop[ p1 ] , mesons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
    } else {
      // I can't think of a time when this would be legitimate
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	printf( "[MESONS] attempt to contract two different source type"
		"propagators thwarted \n" ) ;
	return FAILURE ;
      }
      // check that the two props have the same origin?
      int mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( prop[ p1 ].origin[ mu ] != prop[ p2 ].origin[ mu ] ) {
	  printf( "[MESONS] contraction of mesons with unequal origins"
		  "%d vs %d ( index %d ) " , prop[ p1 ].origin[ mu ] ,
		  prop[ p2 ].origin[ mu ] , mu ) ;
	  return FAILURE ;
	}
      }
      // otherwise we plough on
      if( mesons_offdiagonal( prop[ p1 ] , prop[ p2 ] ,
			      mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
    // loop on measurements
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
