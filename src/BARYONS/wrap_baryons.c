/**
   @file wrap_baryons.c
   @brief baryon contraction wrapper
 */

#include "common.h"

#include "baryons_uuu.h"      // flavour degenerate
#include "baryons_uud.h"      // 2 the same, 1 different
#include "baryons_uds.h"      // 3 different quarks in contraction
#include "GLU_timer.h"        // print_time()
#include "read_propheader.h"  // for read_propheader()

// placeholder for when we do the baryons
int
contract_baryons( struct propagator *prop ,
		  const struct baryon_info *baryons ,
		  const struct cut_info CUTINFO ,
		  const size_t nbaryons )
{
  fprintf( stdout , "\n[BARYONS] performing %zu contraction(s)\n" , nbaryons ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nbaryons ; measurements++ ) {
    // to make it more legible
    const size_t p1 = baryons[ measurements ].map[0] ;
    const size_t p2 = baryons[ measurements ].map[1] ;
    const size_t p3 = baryons[ measurements ].map[2] ;

    // check origins are the same and plaquettes are the same
    if( sanity_check_props( prop , baryons[ measurements ].map ,
			    3 , "[BARYONS]" ) == FAILURE ) {
      return FAILURE ;
    }
    
    // big logic block for contracting the right pieces
    if( p1 == p2 && p2 == p3 ) {
      // only have flavour diagonal option at the moment
      if( baryons_diagonal( prop[ p1 ] , CUTINFO , 
			    baryons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      // rewind file and read header again
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      // two props are the same S3 ( S1 Cgmu S1 Cgmu )
    } else if( ( p1 == p2 && p2 != p3 ) ) {
      if( baryons_2fdiagonal( prop[ p1 ] , prop[ p3 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
      // two props are the same S1 ( S1 Cgmu S2 )
    } else if( p1 == p3 && p3 != p2 ) {
      if( baryons_2fdiagonal( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
      // two props are the same S1 ( S2 Cgmu S2 Cgmu )
    } else if( p2 == p3 && p1 != p2 ) {
      if( baryons_2fdiagonal( prop[ p2 ] , prop[ p1 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      // otherwise we resort to the 3-component baryon
    } else {
      if( baryons_3fdiagonal( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ,
			      CUTINFO , baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
    }
    // tell us how long it took
    print_time( ) ;
  }
  return SUCCESS ;
}
