/**
   @file wrap_baryons.c
   @brief baryon contraction wrapper
 */

#include "common.h"

#include "baryons.h"

// placeholder for when we do the baryons
int
contract_baryons( struct propagator *prop ,
		  const struct meson_info *baryons ,
		  const int nbaryons )
{
  printf( "\n[BARYONS] performing %d contraction(s) \n" , nbaryons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nbaryons ; measurements++ ) {
    // to make it more legible
    const int p1 = baryons[ measurements ].map[0] ;
    const int p2 = baryons[ measurements ].map[1] ;

    // only have flavour diagonal option at the moment
    if( baryons_diagonal( prop[ p1 ] , baryons[ measurements ].outfile ) 
	== FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
