/**
   @file wrap_baryons.c
   @brief baryon contraction wrapper
 */

#include "common.h"

#include "baryons_uuu.h"      // flavour degenerate
#include "baryons_uud.h"      // 2 the same, 1 different
#include "baryons_uds.h"      // 3 different quarks in contraction
#include "GLU_timer.h"        // print_time()
#include "read_propheader.h"  // for read_propheader()

// make sure we have the same source origins
static int
check_origins( struct propagator p1 ,
	       struct propagator p2 ,
	       struct propagator p3 )
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( ( p1.origin[ mu ] != p2.origin[ mu ] ) || 
	( p1.origin[ mu ] != p3.origin[ mu ] ) ||
	( p2.origin[ mu ] != p3.origin[ mu ] ) ) {
      fprintf( stderr , "[BARYONS] contraction of baryons with unequal origins"
	       "%zu vs %zu vs. %zu ( index %zu )\n" , 
	       p1.origin[ mu ] , p2.origin[ mu ] , p3.origin[ mu ] , mu ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// placeholder for when we do the baryons
int
contract_baryons( struct propagator *prop ,
		  const struct baryon_info *baryons ,
		  const struct cut_info CUTINFO ,
		  const size_t nbaryons )
{
  fprintf( stdout , "\n[BARYONS] performing %zu contraction(s)\n" , nbaryons ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nbaryons ; measurements++ ) {
    // to make it more legible
    const size_t p1 = baryons[ measurements ].map[0] ;
    const size_t p2 = baryons[ measurements ].map[1] ;
    const size_t p3 = baryons[ measurements ].map[2] ;

    // big logic block for contracting the right pieces
    if( p1 == p2 && p2 == p3 ) {
      // only have flavour diagonal option at the moment
      if( baryons_diagonal( prop[ p1 ] , CUTINFO , 
			    baryons[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      // rewind file and read header again
      rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      // two props are the same S3 ( S1 Cgmu S1 Cgmu )
    } else if( ( p1 == p2 && p2 != p3 ) ) {
      if( prop[ p1 ].source != prop[ p3 ].source ) {
	fprintf( stderr , "[BARYONS] Caught unequal sources contraction\n" ) ;
	return FAILURE ;
      }
      if( check_origins( prop[ p1 ] , prop[ p1 ] , prop[ p3 ] ) == FAILURE ) {
	return FAILURE ;
      }
      if( baryons_2fdiagonal( prop[ p1 ] , prop[ p3 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
      // two props are the same S2 ( S1 Cgmu S1 Cgmu )
    } else if( p1 == p3 && p3 != p2 ) {
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	fprintf( stderr , "[BARYONS] Caught unequal sources contraction\n" ) ;
	return FAILURE ;
      }
      if( check_origins( prop[ p1 ] , prop[ p1 ] , prop[ p2 ] ) == FAILURE ) {
	return FAILURE ;
      }
      if( baryons_2fdiagonal( prop[ p1 ] , prop[ p2 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
      // two props are the same S1 ( S2 Cgmu S2 Cgmu )
    } else if( p2 == p3 && p1 != p2 ) {
      if( prop[ p2 ].source != prop[ p1 ].source ) {
	fprintf( stderr , "[BARYONS] Caught unequal sources contraction\n" ) ;
	return FAILURE ;
      }
      if( check_origins( prop[ p1 ] , prop[ p1 ] , prop[ p2 ] ) == FAILURE ) {
	return FAILURE ;
      }
      if( baryons_2fdiagonal( prop[ p2 ] , prop[ p1 ] , CUTINFO , 
			      baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
      rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      // otherwise we resort to the 3-component baryon
    } else {
      if( prop[ p1 ].source != prop[ p2 ].source ||
	  prop[ p1 ].source != prop[ p3 ].source ||
	  prop[ p2 ].source != prop[ p3 ].source ) {
	fprintf( stderr , "[BARYONS] Caught unequal sources contraction\n" ) ;
	return FAILURE ;
      }
      if( check_origins( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ) == FAILURE ) {
	return FAILURE ;
      }
      if( baryons_3fdiagonal( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ,
			      CUTINFO , baryons[ measurements ].outfile ) 
	  == FAILURE ) {
	return FAILURE ;
      }
      rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
      rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
    }
    // tell us how long it took
    print_time( ) ;
  }
  return SUCCESS ;
}
