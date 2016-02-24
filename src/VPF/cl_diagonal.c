/**
   @file cl_diagonal.c
   @brief conserved-local Wilson currents

   At the moment I treat the NREL like the chiral which is
   probably wrong, but we cross that bridge when we get to it
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "currents.h"          // conserved-local
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time()
#include "io.h"                // read_prop
#include "matrix_ops.h"        // constant_mul_gauge
#include "momspace_PImunu.h"   // momentum space VPF
#include "read_propheader.h"   // reread the header
#include "tmoments_PImunu.h"   // time moments too
#include "setup.h"             // general setup
#include "spinor_ops.h"        // spinor_minus
#include "WardIdentity.h"      // Config space WI

// number of propagators
#define Nprops (2)

// compute the conserved local for a correlator
int
cl_diagonal( struct propagator prop1 ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile )
{
  // counters
  const size_t stride1 = NS ;
  const size_t stride2 = NS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // vector gamma map
  const size_t VGMAP[ ND ] = { GAMMA_X , GAMMA_Y , GAMMA_Z , GAMMA_T } ;

  // need to look these up
  const size_t AGMAP[ ND ] = { 6 , 7 , 8 , 9 } ;

  // allocations for the PI-data
  struct PIdata *DATA_AA = NULL ; 
  struct PIdata *DATA_VV = NULL ; 

  // loop counters
  size_t x , t ;
  
  // error code
  int error_code = SUCCESS ;

  // initialise our measurement struct
  const struct propagator prop[ Nprops ] = { prop1 , prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = calloc( LVOLUME , sizeof( struct PIdata ) ) ;
  DATA_VV = calloc( LVOLUME , sizeof( struct PIdata ) ) ;

  // Read first timeslice and the one above it
  if( read_prop( prop1 , M.S[0] ) == FAILURE ||
      read_prop( prop1 , M.S[1] ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // copy for the final timeslice
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    equate_spinor_minus( &M.Sf[1][x] , &M.S[0][x] ) ;
  }

  // NB :
  // this timeslice is in M.S[0]
  // upper timeslice is in M.S[1]
  // upper+1 timeslice is in M.Sf[0]
  // final timeslice is in M.Sf[1]

  for( t = 0 ; t < ( LT-1 ) ; t++ ) {

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ ND-1 ] ) % LT ;
    
    // parallel loop with an error flag
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT-2 ) ) {
	  if( read_prop( prop1 , M.Sf[0] ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      #pragma omp for private(x) schedule(dynamic)
      for( x = 0 ; x < LCU ; x++ ) {
	// do the conserved-local contractions
	contract_conserved_local_site( DATA_AA , DATA_VV , 
				       lat , 
				       M.S[0] , M.S[1] , 
				       M.S[0] , M.S[1] ,
				       M.GAMMAS , AGMAP , VGMAP , x , 
				       tshifted ) ;
      }
    }

    // leave if something went bad
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &M.S[0][x] , &M.S[1][x] , sizeof( struct spinor ) ) ;
      memcpy( &M.S[1][x] , &M.Sf[0][x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] cl-flavour diagonal done %.f %%", 
	   (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }

  // and contract the final timeslice
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    contract_conserved_local_site( DATA_AA , DATA_VV , 
				   lat , 
				   M.S[0] , M.Sf[1] , 
				   M.S[0] , M.Sf[1] , 
				   M.GAMMAS , AGMAP , VGMAP , x , 
				   ( t + LT - prop1.origin[ ND-1 ] ) % LT  ) ;
  }
  printf("\r[VPF] cl-flavour diagonal done 100%% \n" ) ; 

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  if( error_code != FAILURE ) {

    // derivatives delta_\mu V_\mu(x)
    WI_configspace_bwd( DATA_VV , lat ) ;

    // compute the time moments
    tmoments( DATA_AA , DATA_VV , outfile , CONSERVED_LOCAL ) ;

    // do all the momspace stuff away from the contractions
    if( prop1.source == POINT ) { 
      momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		       CONSERVED_LOCAL ) ;
    }
  }
  print_time( ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return error_code ;
}

#undef Nprops
