/**
   @file wrap_VPF.c
   @brief wrapper for calling the VPF contractions
 */
#include "common.h"

#include "cl_diagonal.h"      // conserved-local Wilson currents
#include "cl_offdiagonal.h"   // conserved-local flavour off diagonal
#include "GLU_timer.h"        // print_time()
#include "ll_diagonal.h"      // local-local currents
#include "ll_offdiagonal.h"   // flavour off diagonal local-local
#include "read_propheader.h"  // read_propheader()

static int 
( *single_callback ) ( struct propagator prop ,
		       const struct site *lat ,
		       const struct cut_info CUTINFO ,
		       const char *outfile ) ;

static void
select_callback_single( const current_type current )
{
  switch( current ) {
  case CONSERVED_LOCAL :
    single_callback = cl_diagonal ;
    break ;
  case LOCAL_LOCAL :
    single_callback = ll_diagonal ;
    break ;
    // have space for conserved-conserved if we choose to do it
  }
  return ;
}

static int 
( *double_callback ) ( struct propagator prop1 ,
		       struct propagator prop2 ,
		       const struct site *lat ,
		       const struct cut_info CUTINFO ,
		       const char *outfile ) ;

static void
select_callback_double( const current_type current )
{
  switch( current ) {
  case CONSERVED_LOCAL :
    double_callback = cl_offdiagonal ;
    break ;
  case LOCAL_LOCAL :
    double_callback = ll_offdiagonal ;
    break ;
  }
  return ;
}

// meson contraction driver
int
contract_VPF( struct propagator *prop ,
	      const struct site *lat ,
	      const struct VPF_info *VPF ,
	      const struct cut_info CUTINFO ,
	      const size_t nVPF )
{
  fprintf( stdout , "\n[VPF] performing %zu contraction(s) \n" , nVPF ) ;
  size_t measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nVPF ; measurements++ ) {
    const size_t p1 = VPF[ measurements ].map[0] ;
    const size_t p2 = VPF[ measurements ].map[1] ;

    if( p1 == p2 ) {
      select_callback_single( VPF[ measurements ].current ) ;
      // and we use the function pointer we have set
      if( single_callback( prop[ p1 ] , lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      // rewind file and read header again
      rewind( prop[p1].file ) ; read_propheader( &prop[p1] ) ;
    } else {
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	fprintf( stderr , "[VPF] thwarted attempt contracting"
		 " different sources \n" ) ;
	return FAILURE ;
      }
      select_callback_double( VPF[ measurements ].current ) ;
      if( double_callback( prop[ p1 ] , prop[ p2 ] , lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
      // rewind file and read header again
      rewind( prop[p1].file ) ; read_propheader( &prop[p1] ) ;
      rewind( prop[p2].file ) ; read_propheader( &prop[p2] ) ;
    }
    // end of measurement loop
    print_time() ;
  }
  // I would consider getting here a success
  return SUCCESS ;
}
