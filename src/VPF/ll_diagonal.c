/**
   @file ll_diagonal.c
   @brief flavour diagonal local_local current contractions
 */

#include "common.h"

#include "basis_conversions.h" // rotate slice
#include "currents.h"          // contractions
#include "gammas.h"            // gamma matrices
#include "io.h"                // read_prop
#include "matrix_ops.h"        // constant_mul_gauge
#include "momspace_PImunu.h"   // momentum space VPF
#include "progress_bar.h"      // progress_bar()
#include "tmoments_PImunu.h"   // time moments calculations
#include "setup.h"             // initialising and stuff

#include "geometry.h"

// number of propagators
#define Nprops (1)

// compute the conserved local for a correlator
int
ll_diagonal( struct propagator prop1 ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile )
{ 
  // vector gamma map
  const size_t VGMAP[ ND ] = { GAMMA_X , GAMMA_Y , GAMMA_Z , GAMMA_T } ;

  // need to look these up
  const size_t AGMAP[ ND ] = { 6 , 7 , 8 , 9 } ;

  // counters
  const size_t stride1 = ND ;
  const size_t stride2 = ND ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // loop counters
  size_t t = 0 , x ;

  // error code
  int error_code = SUCCESS ;

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = NULL , *DATA_VV = NULL ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[VPF] Failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // initially read a timeslice
#pragma omp parallel
  {
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;
  }
  if( error_code == FAILURE ) {
    goto memfree ;
  }

  // loop the timeslices
  for( t = 0 ; t < LT ; t++ ) {

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ ND-1 ] ) % LT ;

    // do the conserved-local contractions
    #pragma omp parallel
    {
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t ) ;
      }
      #pragma omp for private(x)
      for( x = 0 ; x < LCU ; x++ ) {
	contract_local_local_site( DATA_AA , DATA_VV , M.S[0] , M.S[0] , 
				   M.GAMMAS , AGMAP , VGMAP , x , 
				   tshifted ) ;
      }
    }

    // to err is human
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over
    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      memcpy( &M.S[0][x] , &M.Sf[0][x] , sizeof( struct spinor ) ) ;
    }

    // status
    progress_bar( t , LT ) ;
  }

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  if( error_code != FAILURE ) {
    // temporal data too
    tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;
    
    // do all the momspace stuff away from the contractions
    //if( prop1.source == POINT ) {
      momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		       LOCAL_LOCAL ) ;
      //}
  }

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
