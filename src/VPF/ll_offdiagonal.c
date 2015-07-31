/**
   @file ll_offdiagonal.c
   @brief flavour off-diagonal local_local current contractions
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
#include "read_propheader.h"   // read_propheader()

// do the local-local for different fermion props
int
ll_offdiagonal( struct propagator prop1 ,
		struct propagator prop2 ,
		const struct site *lat ,
		const struct cut_info CUTINFO ,
		const char *outfile )
{
  // vector gamma map
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

  // need to look these up
  const int AGMAP[ ND ] = { GAMMA_5 + 1 , GAMMA_5 + 2 , GAMMA_5 + 3 , GAMMA_5 + 4 } ;

  // and our spinor
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = NULL , *DATA_VV = NULL ;

  // allocate spinors 
  if( corr_malloc( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop1.basis ) == FAILURE ) {
    goto free_failure ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // initially read a timeslice
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto free_failure ;
  }

  // loop the timeslices
  int t , x ;
  for( t = 0 ; t < LT ; t++ ) {

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ ND-1 ] ) % LT ;

    // do the conserved-local contractions
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      if( t < LT-1 ) {
      #pragma omp master
	{
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      #pragma omp single nowait
	{
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      // loop spatial volume
      #pragma omp for private(x)
      for( x = 0 ; x < LCU ; x++ ) {
	contract_local_local_site( DATA_AA , DATA_VV , S1 , S2 , 
				   GAMMAS , AGMAP , VGMAP , x , 
				   tshifted ) ;
      }
    }

    // to err is human
    if( error_flag == FAILURE ) {
      goto free_failure ;
    }

    // copy over
    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      memcpy( &S1[x] , &S1f[x] , sizeof( struct spinor ) ) ;
      memcpy( &S2[x] , &S2f[x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] ll-flavour diagonal done %.f %%", 
	   (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;
  free( S2 ) ; free( S2f ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // temporal data too
  tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  if( prop1.source == POINT ) {
    momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		     LOCAL_LOCAL ) ;
  }

  print_time( ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  return SUCCESS ;

 free_failure :

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;
  free( S2 ) ; free( S2f ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return FAILURE ;
}
