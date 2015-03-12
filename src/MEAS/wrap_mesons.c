/**
   @file wrap_mesons.c
   @brief meson contractions wrappers
 */

#include "common.h"

#include "brutal_mesons.h"   // brutal mesons
#include "conserved_local.h" // conserved-local contractions
#include "mesons.h"          // meson contractions
#include "wall_mesons.h"     // wall mesons

static int 
( *single_callback ) ( struct propagator prop ,
		       const char *outfile ) ;

static void
select_callback_single( const sourcetype source )
{
  switch( source ) {
  case WALL :
    single_callback = wall_mesons ;
    break ;
  case POINT :
    single_callback = single_mesons ;
    break ;
  }
  return ;
}

static int 
( *double_callback ) ( struct propagator prop1 ,
		       struct propagator prop2 ,
		       const char *outfile ) ;

static void
select_callback_double( const sourcetype source )
{
  switch( source ) {
  case WALL :
    double_callback = wall_double_mesons ;
    break ;
  case POINT :
    double_callback = double_mesons ;
    break ;
  }
  return ;
}

// meson contraction driver
int
contract_mesons( struct propagator *prop ,
		 const struct meson_info *mesons ,
		 const int nmesons )
{
  printf( "\n[MESONS] performing %d contraction(s) \n" , nmesons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nmesons ; measurements++ ) {
    // to make it more legible
    const int p1 = mesons[ measurements ].map[0] ;
    const int p2 = mesons[ measurements ].map[1] ;

    if( p1 == p2 ) {
      select_callback_single( prop[ p1 ].source ) ;
      // and we use the function pointer we have set
      if( single_callback( prop[ p1 ] ,
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    } else {
      // I can't think of a time when this would be legitimate
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	printf( "[MESONS] attempt to contract two different source type propagators thwarted \n" ) ;
	return FAILURE ;
      }
      select_callback_double( prop[ p1 ].source ) ;
      if( double_callback( prop[ p1 ] ,
			   prop[ p2 ] ,
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
    // loop on measurements
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
