/**
   @file wrap_mesons.c
   @brief meson contractions wrappers
 */
#include "common.h"

#include "GLU_timer.h"       // print_time()
#include "mesons.h"          // flavour diagonal meson contractions
#include "mesons_offdiag.h"  // flavour off-diagonal meson contractions
#include "read_propheader.h" // (re)read the propagator header

// meson contraction driver
int
contract_mesons( struct propagator *prop ,
		 const struct meson_info *mesons ,
		 const struct cut_info CUTINFO ,
		 const size_t nmesons )
{
  fprintf( stdout , "\n[MESONS] performing %zu contraction(s) \n" , nmesons ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nmesons ; measurements++ ) {
    // to make it more legible
    const size_t p1 = mesons[ measurements ].map[0] ;
    const size_t p2 = mesons[ measurements ].map[1] ;

    // check origins are the same and plaquettes are the same
    if( sanity_check_props( prop , mesons[ measurements ].map ,
			    2 , "[MESONS]" ) == FAILURE ) {
      return FAILURE ;
    }

    // flavour diagonal meson
    if( p1 == p2 ) {
      if( mesons_diagonal( prop[ p1 ] , CUTINFO , 
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
    } else {
      // otherwise we plough on
      if( mesons_offdiagonal( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			      mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
    }
    // loop on measurements
    print_time( ) ;
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
