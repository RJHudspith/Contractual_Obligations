/**
   @file wrap_baryons.c
   @brief baryon contraction wrapper
 */

#include "common.h"

#include "baryons.h"   // flavour degenerate
#include "baryons2.h"  // 2 the same, 1 different
#include "baryons3.h"  // 3 different quarks in contraction

// placeholder for when we do the baryons
int
contract_baryons( struct propagator *prop ,
		  const struct baryon_info *baryons ,
		  const int nbaryons )
{
  printf( "\n[BARYONS] performing %d contraction(s) \n" , nbaryons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nbaryons ; measurements++ ) {
    // to make it more legible
    const int p1 = baryons[ measurements ].map[0] ;
    const int p2 = baryons[ measurements ].map[1] ;
    const int p3 = baryons[ measurements ].map[2] ;

    // big logic block for contracting the right pieces
    if( p1 == p2 && p2 == p3 ) {
      // only have flavour diagonal option at the moment
      if( baryons_diagonal( prop[ p1 ] , baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      // two props are the same S3 ( S1 Cgmu S1 Cgmu )
    } else if( ( p1 == p2 && p2 != p3 ) ) {
      if( prop[ p1 ].source != prop[ p3 ].source ) {
	printf( "[BARYONS] Caught unequal sources contraction \n" ) ;
	return FAILURE ;
      }
      if( baryons_2diagonal( prop[ p1 ] , prop[ p3 ] , 
			     baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      // two props are the same S2 ( S1 Cgmu S1 Cgmu )
    } else if( p1 == p3 && p3 != p2 ) {
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	printf( "[BARYONS] Caught unequal sources contraction \n" ) ;
	return FAILURE ;
      }
      if( baryons_2diagonal( prop[ p1 ] , prop[ p2 ] , 
			     baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      // two props are the same S1 ( S2 Cgmu S2 Cgmu )
    } else if( p2 == p3 && p1 != p2 ) {
      if( prop[ p2 ].source != prop[ p1 ].source ) {
	printf( "[BARYONS] Caught unequal sources contraction \n" ) ;
	return FAILURE ;
      }
      if( baryons_2diagonal( prop[ p2 ] , prop[ p1 ] , 
			     baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      // otherwise we resort to the 3-component baryon
    } else {
      if( prop[ p1 ].source != prop[ p2 ].source ||
	  prop[ p1 ].source != prop[ p3 ].source ||
	  prop[ p2 ].source != prop[ p3 ].source ) {
	printf( "[BARYONS] Caught unequal sources contraction \n" ) ;
	return FAILURE ;
      }
      if( baryons_3fdiagonal( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ,
			     baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
    }
  }
  return SUCCESS ;
}
