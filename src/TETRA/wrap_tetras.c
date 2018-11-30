/**
   @file wrap_tetras.c
   @brief tetraquark contraction wrapper
 */

#include "common.h"

#include "GLU_timer.h"        // print_time()
#include "HAL_su2.h"
#include "su2_dibaryon.h"     // su2_dibaryon()
#include "su2_rhoeta.h"       // su2_rhoeta()
#include "tetra_udbb.h"       // light flavour degenerate
#include "tetra_usbb.h"       // light flavour agnostic heavy degen
#include "tetra_udcb.h"       // light flavour degenerate heavy not
#include "tetra_uscb.h"       // all non-degenerate
#include "read_propheader.h"  // for read_propheader()

//#define SU2_RHOETA
#define HAL_RHORHO

#if NC == 3

// first case :: udbb type there are only two props
static int
contract_udbb( struct propagator *prop ,
	       const struct cut_info CUTINFO ,
	       const char *outfile ,
	       const size_t p1 ,
	       const size_t p2 )
{
  fprintf( stdout , "[TETRA] contracting a udbb type tetra\n" ) ;
  if( tetraquark_udbb( prop[ p1 ] , prop[ p2 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
  return SUCCESS ;
}

// second case :: lights are not degenerate
static int
contract_usbb( struct propagator *prop ,
	       const struct cut_info CUTINFO ,
	       const char *outfile ,
	       const size_t p1 ,
	       const size_t p2 ,
	       const size_t p3 )
{
  fprintf( stdout , "[TETRA] contracting a usbb type tetra\n" ) ;
  if( tetraquark_usbb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
  return SUCCESS ;
}

// third case :: heavys are not degenerate
static int
contract_udcb( struct propagator *prop ,
	       const struct cut_info CUTINFO ,
	       const char *outfile ,
	       const size_t p1 ,
	       const size_t p2 ,
	       const size_t p3 )
{
  fprintf( stdout , "[TETRA] contracting a udcb type tetra\n" ) ;
  if( tetraquark_udcb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
  return SUCCESS ;
}

// fourth case :: all are not degenerate
static int
contract_uscb( struct propagator *prop ,
	       const struct cut_info CUTINFO ,
	       const char *outfile , 
	       const size_t p1 ,
	       const size_t p2 ,
	       const size_t p3 ,
	       const size_t p4 )
{
  fprintf( stdout , "[TETRA] contracting a uscb type tetra\n" ) ;
  if( tetraquark_uscb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , prop[ p4 ] , 
		       CUTINFO , outfile )
      == FAILURE ) {
    return FAILURE ;
  }
  if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
  if( reread_propheaders( &prop[ p4 ] ) == FAILURE ) { return FAILURE ; }
  return SUCCESS ;
}

#endif

// tetraquark calculator, prop3 should be the heavy one
int
contract_tetras( struct propagator *prop ,
		 const struct tetra_info *tetras ,
		 const struct cut_info CUTINFO ,
		 const size_t ntetras )
{
  fprintf( stdout , "\n[TETRA] performing %zu contraction(s) \n" , ntetras ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < ntetras ; measurements++ ) {
    // to make it more legible
    const size_t p1 = tetras[ measurements ].map[0] ;
    const size_t p2 = tetras[ measurements ].map[1] ;
    const size_t p3 = tetras[ measurements ].map[2] ;
    const size_t p4 = tetras[ measurements ].map[3] ;

    // check origins are the same and plaquettes are the same
    if( sanity_check_props( prop , tetras[ measurements ].map ,
			    4 , "[TETRAS]" ) == FAILURE ) {
      return FAILURE ;
    }

    #if NC == 2

    fprintf( stdout , "[TETRA] su2 dibaryon contractions\n" ) ;

    if( p1 == p2 && p2 == p3 && p3 == p4 ) {

      #ifdef SU2_RHOETA
      if( su2_rhoeta( prop[ p1 ] , CUTINFO , tetras[ measurements ].outfile )
	  == FAILURE ) {
	return FAILURE ;
      }
      #elif (defined HAL_RHORHO)
      if( HAL_su2( prop[ p1 ] , CUTINFO , tetras[ measurements ].outfile )
	  == FAILURE ) {
	return FAILURE ;
      }      
      #else
      // su2 dibaryon code
      if( su2_dibaryon( prop[ p1 ] , CUTINFO , tetras[ measurements ].outfile )
	  == FAILURE ) {
	return FAILURE ;
      }
      #endif
      if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
    } else {
      fprintf( stderr , "[TETRA] non-similar case not supported\n" ) ;
      return FAILURE ;
    }

    #elif NC == 3

    // support for degenerate light content
    if( p1 == p2 ) {
      // sanity trap
      if( p2 == p3 || p2 == p4 || p1 == p3 || p1 == p4 ) {
	fprintf( stderr , "Tetraquark ( %zu %zu %zu %zu ) not supported \n" , 
		 p1 , p2 , p3 , p4 ) ;
	return FAILURE ;
      }
      // degenerate heavy quarks
      if( p3 == p4 ) {
	// degenerate light and heavies
	if( contract_udbb( prop , CUTINFO , tetras[ measurements ].outfile ,
			   p1 , p3 ) == FAILURE ) {
	  return FAILURE ;
	}
      } else {
	// non-degenerate heavies
	if( contract_udcb( prop , CUTINFO , tetras[ measurements ].outfile ,
			   p1 , p3 , p4 ) == FAILURE ) {
	  return FAILURE ;
	}
      }
    } else if( p3 == p4 ) {
      // sanity trap
      if( p2 == p3 || p2 == p4 || p1 == p3 || p1 == p4 ) {
	printf( "Tetraquark ( %zu %zu %zu %zu ) not supported \n" , 
		p1 , p2 , p3 , p4 ) ;
      }
      // degenerate light quarks already covered above so only have
      // usbb case
      if( contract_usbb( prop , CUTINFO , tetras[ measurements ].outfile ,
			 p1 , p2 , p3 ) == FAILURE ) { 
	return FAILURE ;
      }
    } else {
      // all non degenerate eats all the other options
      if( contract_uscb( prop , CUTINFO , tetras[ measurements ].outfile ,
			 p1 , p2 , p3 , p4 ) == FAILURE ) {
	return FAILURE ;
      }
      //
    }

    #else
    frpintf( stderr , "[TETRA] tetra contractions not supported for NC = %d\n" ,
	     NC ) ;
    return FAILURE ;
    #endif
    
    // tell us how long it all took
    print_time( ) ;
  }
  return SUCCESS ;
}
