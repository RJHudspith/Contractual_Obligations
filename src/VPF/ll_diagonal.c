/**
   @file ll_diagonal.c
   @brief flavour diagonal local_local current contractions
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

// compute the conserved local for a correlator
int
ll_diagonal( struct propagator prop ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile )
{
  // vector gamma map
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

  // need to look these up
  const int AGMAP[ ND ] = { GAMMA_5 + 1 , GAMMA_5 + 2 , GAMMA_5 + 3 , GAMMA_5 + 4 } ;

  // and our spinor
  struct spinor *S1 = NULL , *S1f = NULL ;

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

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    goto free_failure ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // initially read a timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto free_failure ;
  }

  // loop the timeslices
  int t , x ;
  for( t = 0 ; t < L0 ; t++ ) {

    // do the conserved-local contractions
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < L0-1 ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      #pragma omp for private(x)
      for( x = 0 ; x < LCU ; x++ ) {
	contract_local_local_site( DATA_AA , DATA_VV , S1 , S1 , 
				   GAMMAS , AGMAP , VGMAP , x , t ) ;
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
    }

    // status
    printf("\r[VPF] ll-flavour diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // temporal data too
  tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  if( prop.source == POINT ) {
    momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		     LOCAL_LOCAL ) ;
  }

  print_time( ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  return SUCCESS ;

 free_failure :

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return FAILURE ;
}
