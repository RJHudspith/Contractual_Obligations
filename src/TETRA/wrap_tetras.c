/**
   @file wrap_tetras.c
   @brief tetraquark contraction wrapper
 */

#include "common.h"

#include "tetraquark.h"   // flavour degenerate

// for origin checking
static int
check_origins( struct propagator p1 ,
	       struct propagator p2 ,
	       struct propagator p3 ,
	       struct propagator p4 )
{
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( p1.origin[mu] != p2.origin[mu] ||
	p1.origin[mu] != p3.origin[mu] ||
	p1.origin[mu] != p4.origin[mu] ||
	p2.origin[mu] != p3.origin[mu] ||
	p2.origin[mu] != p4.origin[mu] ||
	p3.origin[mu] != p4.origin[mu] ) {
      printf( "[TETRA] mismatched origins %d %d %d %d (index %d)\n" ,
	      p1.origin[mu] , p2.origin[mu] ,
	      p3.origin[mu] , p4.origin[mu] , mu ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

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

    if( p3 == p4 ) {
      if( tetraquark( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , CUTINFO , 
		      tetras[ measurements ].outfile
		      ) == FAILURE ) {
	    return FAILURE ;
	  }
    } else {
      if( check_origins( prop[ p1 ] , prop[ p2 ] , 
			 prop[ p3 ] , prop[ p4 ] ) == FAILURE ) {
	return FAILURE ;
      }
      printf( "[TETRA] non-degenerate not supported\n" ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
