/**
   @file wrap_diquark.c
   @brief wrapper for the diquark contractions
 */
#include "common.h"

// diquark contractions
#include "diquark.h"
#include "diquark_degen.h"

// contract the diquarks
int
contract_diquarks( struct propagator *prop ,
		   const struct meson_info *diquarks ,
		   const struct cut_info CUTINFO ,
		   const size_t ndiquarks )
{
  printf( "\n[DIQUARK] performing %zu contraction(s) \n" , ndiquarks ) ;
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
    } else {
      if( diquark( prop[ p1 ] , prop[ p1 ] , CUTINFO , 
		   diquarks[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
  }
  return SUCCESS ;
}
