/**
   @file ll_offdiagonal.c
   @brief flavour off-diagonal local_local current contractions

   TODO :: nrel rotation?
 */
#include "common.h"

#include "basis_conversions.h" // rotate slice
#include "currents.h"          // contractions
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time()
#include "io.h"                // read_prop
#include "matrix_ops.h"        // constant_mul_gauge
#include "momspace_PImunu.h"   // momentum space VPF
#include "tmoments_PImunu.h"   // time moments calculations
#include "setup.h"             // init_measurements()

// number of propagators
#define Nprops (2)

// do the local-local for different fermion props
int
ll_offdiagonal( struct propagator prop1 ,
		struct propagator prop2 ,
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
  size_t x , t ;

  // error code
  int error_code = SUCCESS ;

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = NULL , *DATA_VV = NULL ;

  // initialise our measurement struct
  const struct propagator prop[ Nprops ] = { prop1 , prop2 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }
 
  // over the whole volume is not as expensive as you may think
  DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // initially read a timeslice
  if( read_prop( prop1 , M.S[0] ) == FAILURE ||
      read_prop( prop2 , M.S[1] ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // loop the timeslices
  for( t = 0 ; t < LT ; t++ ) {

    // rotate if we must
    rotate_offdiag( M.S , prop , Nprops ) ;

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ ND-1 ] ) % LT ;

    // do the conserved-local contractions
    #pragma omp parallel
    {
      if( t < LT-1 ) {
      #pragma omp master
	{
	  if( read_prop( prop1 , M.Sf[0] ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      #pragma omp single nowait
	{
	  if( read_prop( prop2 , M.Sf[1] ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      // loop spatial volume
      #pragma omp for private(x)
      for( x = 0 ; x < LCU ; x++ ) {
	contract_local_local_site( DATA_AA , DATA_VV , 
				   M.S[0] , M.S[1] , 
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
      memcpy( &M.S[1][x] , &M.Sf[1][x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] ll-flavour diagonal done %.f %%", 
	   (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // deallocs
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  if( error_code != FAILURE ) {
    // temporal data too
    tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;
    
    // do all the momspace stuff away from the contractions
    if( prop1.source == POINT ) {
      momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		       LOCAL_LOCAL ) ;
    }
  }
  print_time( ) ;

  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return error_code ;
}

// clean up the number of props
#undef Nprops
