/**
   @file wrap_tetra.c
   @brief tetraquark contraction wrapper
 */

#include "common.h"

#include "tetraquark.h"   // flavour degenerate

// placeholder for when we do the baryons
int
contract_tetras( struct propagator *prop ,
		 const struct tetra_info *tetras ,
		 const struct cut_info CUTINFO ,
		 const int ntetras )
{
  printf( "\n[TETRA] performing %d contraction(s) \n" , ntetras ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < ntetras ; measurements++ ) {
    // to make it more legible
    const int p1 = tetras[ measurements ].map[0] ;
    const int p2 = tetras[ measurements ].map[1] ;
    const int p3 = tetras[ measurements ].map[2] ;
    const int p4 = tetras[ measurements ].map[3] ;

    if( p1 == p2 && p2 == p3 && p3 == p4 ) {
      if( tetraquark_diagonal( prop[ p1 ] , CUTINFO , 
			       tetras[ measurements ].outfile
			       ) == FAILURE ) {
	    return FAILURE ;
	  }
    } else {
      printf( "[TETRA] non-degenerate not supported\n" ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
