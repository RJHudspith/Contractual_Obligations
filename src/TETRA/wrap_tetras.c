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

// for origin checking
static int
check_origins( struct propagator p1 ,
	       struct propagator p2 ,
	       struct propagator p3 , 
	       struct propagator p4 )
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( p1.origin[mu] != p2.origin[mu] ||
	p1.origin[mu] != p3.origin[mu] ||
	p1.origin[mu] != p4.origin[mu] ||
	p2.origin[mu] != p3.origin[mu] || 
	p2.origin[mu] != p4.origin[mu] ||
	p3.origin[mu] != p4.origin[mu] ) {
      fprintf( stderr , "[TETRA] mismatched origins "
	       "%zu %zu %zu %zu (index %zu)\n" ,
	      p1.origin[mu] , p2.origin[mu] , p3.origin[mu] ,
	      p4.origin[mu] , mu ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// first case :: udbb type there are only two props
static int
contract_udbb( struct propagator *prop ,
	       const struct cut_info CUTINFO ,
	       const char *outfile ,
	       const size_t p1 ,
	       const size_t p2 )
{
  fprintf( stdout , "[TETRA] contracting a udbb type tetra\n" ) ;
  if( check_origins( prop[ p1 ] , prop[ p1 ] , 
		     prop[ p2 ] , prop[ p2 ] ) == FAILURE ) {
    return FAILURE ;
  }
  if( tetraquark_udbb( prop[ p1 ] , prop[ p2 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( prop[ p1 ].basis != NREL_CORR ) {
    rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
  }
  if( prop[ p2 ].basis != NREL_CORR ) {
    rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
  }
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
  if( check_origins( prop[ p1 ] , prop[ p2 ] , 
		     prop[ p3 ] , prop[ p3 ] ) == FAILURE ) {
    return FAILURE ;
  }
  if( tetraquark_usbb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( prop[ p1 ].basis != NREL_CORR ) {
    rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
  }
  if( prop[ p2 ].basis != NREL_CORR ) {
    rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
  }
  if( prop[ p3 ].basis != NREL_CORR ) {
    rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
  }
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
  if( check_origins( prop[ p1 ] , prop[ p1 ] , 
		     prop[ p2 ] , prop[ p3 ] ) == FAILURE ) {
    return FAILURE ;
  }
  if( tetraquark_udcb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
		       CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  if( prop[ p1 ].basis != NREL_CORR ) {
    rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
  }
  if( prop[ p2 ].basis != NREL_CORR ) {
    rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
  }
  if( prop[ p3 ].basis != NREL_CORR ) {
    rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
  }
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
  if( check_origins( prop[ p1 ] , prop[ p2 ] , 
		     prop[ p3 ] , prop[ p4 ] ) == FAILURE ) {
    return FAILURE ;
  }
  if( tetraquark_uscb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , prop[ p4 ] , 
		       CUTINFO , outfile )
      == FAILURE ) {
    return FAILURE ;
  }
  if( prop[ p1 ].basis != NREL_CORR ) {
    rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
  }
  if( prop[ p2 ].basis != NREL_CORR ) {
    rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
  }
  if( prop[ p3 ].basis != NREL_CORR ) {
    rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
  }
  if( prop[ p4 ].basis != NREL_CORR ) {
    rewind( prop[ p4 ].file ) ; read_propheader( &prop[ p4 ] ) ;
  }
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
      if( prop[ p1 ].basis != NREL_CORR ) {
	rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
      }
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
