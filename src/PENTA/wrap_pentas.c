/**
   @file wrap_pentas.c
   @brief pentaquark contraction wrapper
 */
#include "common.h"

#include "GLU_timer.h"        // print_time()
#include "penta_udusb.h"      // light flavour degenerate
#include "read_propheader.h"  // for read_propheader()

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
	p2.origin[mu] != p3.origin[mu]  ) {
      fprintf( stderr , "[PENTA] mismatched origins "
	       "%zu %zu %zu (index %zu)\n" ,
	      p1.origin[mu] , p2.origin[mu] , p3.origin[mu] , mu ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// only case three propagators
static int
contract_udusb( struct propagator *prop ,
		const struct cut_info CUTINFO ,
		const char *outfile ,
		const size_t p1 ,
		const size_t p2 ,
		const size_t p3 )
{
  fprintf( stdout , "[PENTA] contracting a udusb type pentaquark\n" ) ;
  if( check_origins( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] ) == FAILURE ) {
    return FAILURE ;
  }
  if( pentaquark_udusb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
			CUTINFO , outfile ) == FAILURE ) {
    return FAILURE ;
  }
  rewind( prop[ p1 ].file ) ; read_propheader( &prop[ p1 ] ) ;
  rewind( prop[ p2 ].file ) ; read_propheader( &prop[ p2 ] ) ;
  rewind( prop[ p3 ].file ) ; read_propheader( &prop[ p3 ] ) ;
  return SUCCESS ;
}

// pentaaquark calculator, prop3 should be the heavy one, prop2 the strange
// and prop1 the light ones
int
contract_pentas( struct propagator *prop ,
		 const struct penta_info *pentas ,
		 const struct cut_info CUTINFO ,
		 const size_t npentas )
{
  printf( "\n[PENTA] performing %zu contraction(s) \n" , npentas ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < npentas ; measurements++ ) {
    // to make it more legible
    const size_t p1 = pentas[ measurements ].map[0] ;
    const size_t p2 = pentas[ measurements ].map[1] ;
    const size_t p3 = pentas[ measurements ].map[2] ;

    // pentaquark
    if( contract_udusb( prop , CUTINFO , pentas[ measurements ].outfile ,
			p1 , p2 , p3 ) == FAILURE ) {
      return FAILURE ;
    }
    
    // tell us how long it all took
    print_time( ) ;
  }
  return SUCCESS ;
}
