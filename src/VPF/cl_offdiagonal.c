/**
   @file cl_offdiagonal.c
   @brief flavour off-diagonal conserved-local Wilson currents
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
#include "spinor_ops.h"        // spinor_minus
#include "WardIdentity.h"      // Config space WI

// compute the conserved local for a correlator
int
cl_offdiagonal( struct propagator prop1 ,
		struct propagator prop2 ,
		const struct site *lat ,
		const struct cut_info CUTINFO ,
		const char *outfile )
{
  // vector gamma map -> x,y,z,t
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

  // need to look these up
  const int AGMAP[ ND ] = { GAMMA_5 + 1 , GAMMA_5 + 2 , GAMMA_5 + 3 , GAMMA_5 + 4 } ;

  // spinor storage
  struct spinor *S1 = NULL , *S1f = NULL , *S1END = NULL , *S1UP = NULL ;
  struct spinor *S2 = NULL , *S2f = NULL , *S2END = NULL , *S2UP = NULL ;

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // finally our result
  struct PIdata *DATA_AA = NULL , *DATA_VV = NULL ;

  // spinors
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S1UP , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S1END , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2UP , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S2END , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // free everything if the read fails
  if( read_prop( prop1 , S1 ) == FAILURE || 
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto free_failure ;
  }

  // if we are doing nonrel-chiral mesons we switch chiral to nrel
  if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
    nrel_rotate_slice( S1 ) ;
  } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
    nrel_rotate_slice( S2 ) ;
  }

  if( read_prop( prop1 , S1UP ) == FAILURE || 
      read_prop( prop2 , S2UP ) == FAILURE ) {
    goto free_failure ;
  }

  // copy for the final timeslice
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    equate_spinor_minus( &S1END[x] , &S1[x] ) ;
    equate_spinor_minus( &S2END[x] , &S2[x] ) ;
  }

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0-1 ; t++ ) {

    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < L0-2 ) {
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	  // rotate as necessary
	  if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
	    nrel_rotate_slice( S1UP ) ;
	  }
	}
      }
      #pragma omp single nowait
      {
	if( t < L0-2 ) {
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	  // rotate as necessary
	  if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
	    nrel_rotate_slice( S2f ) ;
	  }
	}
      }
      #pragma omp for private(x) schedule(dynamic)
      for( x = 0 ; x < LCU ; x++ ) {
	// do the conserved-local contractions
	contract_conserved_local_site( DATA_AA , DATA_VV , 
				       lat , S1 , S1UP , S2 , S2UP ,
				       GAMMAS , AGMAP , VGMAP , x , t ) ;
      }
    }

    // to err is human
    if( error_flag == FAILURE ) {
      goto free_failure ;
    }

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
      memcpy( &S2[x] , &S2UP[x] , sizeof( struct spinor ) ) ;
      memcpy( &S1UP[x] , &S1f[x] , sizeof( struct spinor ) ) ;
      memcpy( &S2UP[x] , &S2f[x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] cl-flavour diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }

  // and contract the final timeslice
#pragma omp for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    contract_conserved_local_site( DATA_AA , DATA_VV , 
				   lat , S1 , S1END , S2 , S2END ,
				   GAMMAS , AGMAP , VGMAP , x , t ) ;
  }

  printf("\r[VPF] cl-flavour off diagonal done 100%% \n" ) ; 

  // derivatives delta_\mu V_\mu(x)
  WI_configspace_bwd( DATA_VV , lat ) ;

  // free our spinors
  free( S1f ) ; free( S1 ) ; free( S1UP ) ; free( S1END ) ;
  free( S2f ) ; free( S2 ) ; free( S2UP ) ; free( S2END ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // time moments are interesting also
  tmoments( DATA_AA , DATA_VV , outfile , CONSERVED_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  if( prop1.source == POINT ) { 
    momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		     CONSERVED_LOCAL ) ;
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

  // free our spinors
  free( S1f ) ; free( S1 ) ; free( S1UP ) ; free( S1END ) ;
  free( S2f ) ; free( S2 ) ; free( S2UP ) ; free( S2END ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return SUCCESS ;
}