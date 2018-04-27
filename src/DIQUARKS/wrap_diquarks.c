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
    if( p1 == p2 ) {
      if( diquark_degen( prop[ p1 ] , CUTINFO , 
			 diquarks[ measurements ].outfile
			 ) == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[p1].file ) ; read_propheader( &prop[p1] ) ;
    } else {
      if( prop[p1].source != prop[p2].source ) {
	fprintf( stderr , "[DIQUARK] unequal source types"
		 "for the two diquarks\n" ) ;
	return FAILURE ;
      }
      if( diquark_offdiag( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			   diquarks[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[p1].file ) ; read_propheader( &prop[p1] ) ;
      rewind( prop[p2].file ) ; read_propheader( &prop[p2] ) ;
    }
    // give us the time
    print_time( ) ;
  }
  return SUCCESS ;
}
