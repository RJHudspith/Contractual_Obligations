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

    // flavour diagonal meson
    if( p1 == p2 ) {
      if( mesons_diagonal( prop[ p1 ] , CUTINFO , 
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      if( prop[ p1 ].basis != NREL_CORR ) {
	rewind( prop[ p1 ].file ) ;
	read_propheader( &prop[ p1 ] ) ;
      }
    } else {
      // check that the two props have the same origin?
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( prop[ p1 ].origin[ mu ] != prop[ p2 ].origin[ mu ] ) {
	  fprintf( stderr , "[MESONS] contraction of mesons with unequal "
		   "origins %zu vs %zu ( index %zu ) " , 
		   prop[ p1 ].origin[ mu ] , prop[ p2 ].origin[ mu ] , mu ) ;
	  return FAILURE ;
	}
      }
      // otherwise we plough on
      if( mesons_offdiagonal( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			      mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      if( prop[p1].basis != NREL_CORR ) {
	rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      }
      if( prop[p2].basis != NREL_CORR ) {
	rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
      }
    }
    // loop on measurements
    print_time( ) ;
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
