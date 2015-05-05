/**
   @file wrap_disprel.c
   @brief dispersion relation wrapper
 */
#include "common.h"

#include "dispersions.h"  // mesonic dispersion relation

// meson dispersion relation computation
int
contract_disprels( struct propagator *prop ,
		   const struct dispersion_info *dispersions ,
		   const struct cut_info CUTINFO ,
		   const int ndispersions )
{
  printf( "\n[DISPREL] performing %d (ND-1) momentum projected contraction(s) \n" , ndispersions ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < ndispersions ; measurements++ ) {
    // to make it more legible
    const int p1 = dispersions[ measurements ].map[0] ;
    const int p2 = dispersions[ measurements ].map[1] ;

    if( p1 == p2 ) {
      if( dispersions_diagonal( prop[ p1 ] , CUTINFO ,
				dispersions[ measurements ].outfile )
	  == FAILURE ) {
	return FAILURE ;
      }
    } else {
      // I can't think of a time when this would be legitimate
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	printf( "[DISPREL] attempt to contract two different source type"
		"propagators thwarted \n" ) ;
	return FAILURE ;
      }
      // check that the two props have the same origin?
      int mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( prop[ p1 ].origin[ mu ] != prop[ p2 ].origin[ mu ] ) {
	  printf( "[DISPREL] contraction of dispersions with unequal origins"
		  "%d vs %d ( index %d ) " , prop[ p1 ].origin[ mu ] ,
		  prop[ p2 ].origin[ mu ] , mu ) ;
	  return FAILURE ;
	}
      }
      // otherwise we plough on
      if( dispersions_offdiagonal( prop[ p1 ] , prop[ p2 ] , CUTINFO ,
				   dispersions[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
    }
    // loop on measurements
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
