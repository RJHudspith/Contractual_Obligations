/**
   @file wrap_VPF.c
   @brief wrapper for calling the VPF contractions
 */

#include "common.h"

#include "cl_diagonal.h"      // conserved-local Wilson currents
#include "cl_offdiagonal.h"   // conserved-local flavour off diagonal
#include "ll_diagonal.h"      // local-local currents
#include "ll_offdiagonal.h"   // flavour off diagonal local-local

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
	      const int nVPF ,
	      const struct cut_info CUTINFO )
{
  printf( "\n[VPF] performing %d contraction(s) \n" , nVPF ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nVPF ; measurements++ ) {
    const int p1 = VPF[ measurements ].map[0] ;
    const int p2 = VPF[ measurements ].map[1] ;

    if( p1 == p2 ) {
      select_callback_single( VPF[ measurements ].current ) ;
      // and we use the function pointer we have set
      if( single_callback( prop[ p1 ] , lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    } else {
      if( prop[ p1 ].source != prop[ p2 ].source ) {
	printf( "[VPF] thwarted attempt contracting different sources \n" ) ;
	return FAILURE ;
      }
      select_callback_double( VPF[ measurements ].current ) ;
      if( double_callback( prop[ p1 ] , prop[ p2 ] , lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
    // end of measurement loop
  }
  // I would consider getting here a success
  return SUCCESS ;
}
