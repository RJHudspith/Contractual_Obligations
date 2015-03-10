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
( *single_callback ) ( FILE *prop1 , const proptype proptype1 ,
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
( *double_callback ) ( FILE *prop1 , const proptype proptype1 ,
		       FILE *prop2 , const proptype proptype2 ,
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
contract_mesons( FILE **fprops ,
		 const struct meson_info *mesons ,
		 const int nmesons )
{
  printf( "\n[MESONS] performing %d contraction(s) \n" , nmesons ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nmesons ; measurements++ ) {
    if( mesons[ measurements ].map[0] == mesons[ measurements ].map[0] ) {
      // logic bomb, cannot contract same prop with different bases
      if( mesons[ measurements ].proptype1 != mesons[ measurements ].proptype2 ) {
	printf( "[MESONS] we cannot contract a meson with itself with two "
		"different bases \n" ) ;
	return FAILURE ;
      }
      select_callback_single( mesons[ measurements ].source ) ;
      // and we use the function pointer we have set
      if( single_callback( fprops[ mesons[ measurements ].map[0] ] , 
			   mesons[ measurements ].proptype1 ,
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    } else {
      select_callback_double( mesons[ measurements ].source ) ;
      if( double_callback( fprops[ mesons[ measurements ].map[0] ] , 
			   mesons[ measurements ].proptype1 ,
			   fprops[ mesons[ measurements ].map[1] ] , 
			   mesons[ measurements ].proptype2 ,
			   mesons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
    // loop on measurements
  }
  // I would consider getting here to be most successful, quite.
  return SUCCESS ;
}
