/**
   @file wrap_tetras.c
   @brief tetraquark contraction wrapper
 */

#include "common.h"

#include "tetra_degen.h"  // light flavour degenerate
#include "tetraquark.h"   // light flavour agnostic

// for origin checking
static int
check_origins( struct propagator p1 ,
	       struct propagator p2 ,
	       struct propagator p3 )
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( p1.origin[mu] != p2.origin[mu] ||
	p1.origin[mu] != p3.origin[mu] ||
	p2.origin[mu] != p3.origin[mu] ) {
      printf( "[TETRA] mismatched origins %d %d %d (index %zu)\n" ,
	      p1.origin[mu] , p2.origin[mu] ,
	      p3.origin[mu] , mu ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// tetraquark calculator, prop3 should be the heavy one
int
contract_tetras( struct propagator *prop ,
		 const struct tetra_info *tetras ,
		 const struct cut_info CUTINFO ,
		 const size_t ntetras )
{
  printf( "\n[TETRA] performing %zu contraction(s) \n" , ntetras ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < ntetras ; measurements++ ) {
    // to make it more legible
    const size_t p1 = tetras[ measurements ].map[0] ;
    const size_t p2 = tetras[ measurements ].map[1] ;
    const size_t p3 = tetras[ measurements ].map[2] ;

    // support for degenerate light content
    if( p1 == p2 && p2 != p3 ) {
      if( tetraquark_degen( prop[ p1 ] , prop[ p3 ] , CUTINFO , 
			    tetras[ measurements ].outfile
			    ) == FAILURE ) {
	return FAILURE ;
      }
    // general p1 != p2 combination
    } else if( p1 != p2 && p2 != p3 ) {
      if( tetraquark( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , CUTINFO , 
		      tetras[ measurements ].outfile
		      ) == FAILURE ) {
	return FAILURE ;
      }
    // otherwise complain 
    } else {
      if( check_origins( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      printf( "[TETRA] non-degenerate not supported\n" ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
