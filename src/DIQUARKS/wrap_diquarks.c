/**
   @file wrap_diquarks.c
   @brief wrapper for the diquark contractions
 */
#include "common.h"

// diquark contractions
#include "diquark.h"
#include "diquark_degen.h"
#include "GLU_timer.h"        // print_time()
#include "read_propheader.h"

// contract the diquarks
int
contract_diquarks( struct propagator *prop ,
		   const struct meson_info *diquarks ,
		   const struct cut_info CUTINFO ,
		   const size_t ndiquarks )
{
  fprintf( stdout , "\n[DIQUARK] performing %zu "
	   "contraction(s)\n" , ndiquarks ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < ndiquarks ; measurements++ ) {
    // to make it more legible
    const size_t p1 = diquarks[ measurements ].map[0] ;
    const size_t p2 = diquarks[ measurements ].map[1] ;

    // check origins are the same and plaquettes are the same
    if( sanity_check_props( prop , diquarks[ measurements ].map ,
			    2 , "[DIQUARKS]" ) == FAILURE ) {
      return FAILURE ;
    }
    if( p1 == p2 ) {
      if( diquark_degen( prop[ p1 ] , CUTINFO , 
			 diquarks[ measurements ].outfile
			 ) == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
    } else {
      if( diquark_offdiag( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			   diquarks[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
      if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
    }
    // give us the time
    print_time( ) ;
  }
  return SUCCESS ;
}
