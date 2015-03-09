/**
   @file wrap_VPF.c
   @brief wrapper for calling the VPF contractions
 */

#include "common.h"

#include "conserved_local.h"  // conserved-local Wilson currents
#include "local_local.h"      // local-local currents

static int 
( *single_callback ) ( FILE *prop1 , 
		       const proptype proptype1 ,
		       const struct site *lat ,
		       const struct cut_info CUTINFO ,
		       const char *outfile ) ;

static void
select_callback_single( const current_type current )
{
  switch( current ) {
  case CONSERVED_LOCAL :
    single_callback = conserved_local ;
    break ;
  case LOCAL_LOCAL :
    single_callback = local_local ;
    break ;
    // have space for conserved-conserved if we choose to do it
  }
  return ;
}

static int 
( *double_callback ) ( FILE *prop1 , const proptype proptype1 ,
		       FILE *prop2 , const proptype proptype2 ,
		       const struct site *lat ,
		       const struct cut_info CUTINFO ,
		       const char *outfile ) ;

static void
select_callback_double( const current_type current )
{
  switch( current ) {
  case CONSERVED_LOCAL :
    double_callback = conserved_local_double ;
    break ;
  case LOCAL_LOCAL :
    double_callback = local_local_double ;
    break ;
  }
  return ;
}

// meson contraction driver
int
contract_VPF( FILE **fprops ,
	      const struct site *lat ,
	      const struct VPF_info *VPF ,
	      const int nVPF ,
	      const struct cut_info CUTINFO )
{
  printf( "\n[VPF] performing %d contractions \n" , nVPF ) ;
  int measurements ;
  // loops measurements and use mesons information to perform contractions
  for( measurements = 0 ; measurements < nVPF ; measurements++ ) {
    if( VPF[ measurements ].map[0] == VPF[ measurements ].map[0] ) {
      select_callback_single( VPF[ measurements ].current ) ;
      // and we use the function pointer we have set
      if( single_callback( fprops[ VPF[ measurements ].map[0] ] , 
			   VPF[ measurements ].proptype1 ,
			   lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    } else {
      select_callback_double( VPF[ measurements ].current ) ;
      if( double_callback( fprops[ VPF[ measurements ].map[0] ] , 
			   VPF[ measurements ].proptype1 ,
			   fprops[ VPF[ measurements ].map[1] ] , 
			   VPF[ measurements ].proptype2 ,
			   lat , CUTINFO ,
			   VPF[ measurements ].outfile ) == FAILURE ) {
	return FAILURE ;
      }
    }
    // end of measurement loop
  }
  // I would consider getting here a success
  return SUCCESS ;
}
