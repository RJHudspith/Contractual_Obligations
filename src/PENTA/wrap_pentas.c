/**
   @file wrap_pentas.c
   @brief pentaquark contraction wrapper
 */
#include "common.h"

#include "GLU_timer.h"        // print_time()
#include "penta_udusb.h"      // light flavour degenerate
#include "read_propheader.h"  // for read_propheader()

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

    // check origins are the same and plaquettes are the same
    if( sanity_check_props( prop , pentas[ measurements ].map ,
			    3 , "[PENTAS]" ) == FAILURE ) {
      return FAILURE ;
    }

    // pentaquark
    fprintf( stdout , "[PENTA] contracting a udusb type pentaquark\n" ) ;
    if( pentaquark_udusb( prop[ p1 ] , prop[ p2 ] , prop[ p3 ] , 
			  CUTINFO , pentas[ measurements ].outfile )
	== FAILURE ) {
      return FAILURE ;
    }
    if( reread_propheaders( &prop[ p1 ] ) == FAILURE ) { return FAILURE ; }
    if( reread_propheaders( &prop[ p2 ] ) == FAILURE ) { return FAILURE ; }
    if( reread_propheaders( &prop[ p3 ] ) == FAILURE ) { return FAILURE ; }
        
    // tell us how long it all took
    print_time( ) ;
  }
  return SUCCESS ;
}
