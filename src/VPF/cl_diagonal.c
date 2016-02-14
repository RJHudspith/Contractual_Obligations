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
#include "spinor_ops.h"        // spinor_minus
#include "WardIdentity.h"      // Config space WI

// compute the conserved local for a correlator
int
cl_diagonal( struct propagator prop ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile )
{
  // vector gamma map
  const size_t VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

  // need to look these up
  const size_t AGMAP[ ND ] = { 6 , 7 , 8 , 9 } ;

  // spinor memory
  struct spinor *S1f = NULL , *S1 = NULL , *S1UP = NULL , *S1END = NULL ;
 
  // allocate gamma matrices
  struct gamma *GAMMAS = NULL ;

  // allocations for the PI-data
  struct PIdata *DATA_AA = NULL ; 
  struct PIdata *DATA_VV = NULL ; 

  // loop counters
  size_t x , t ;
  
  // error code
  int error_code = SUCCESS ;

  // and our spinor(s)
  if( corr_malloc( (void**)&S1f   , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1    , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1UP  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1END , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // over the whole volume is not as expensive as you may think
  DATA_AA = calloc( LVOLUME , sizeof( struct PIdata ) ) ;
  DATA_VV = calloc( LVOLUME , sizeof( struct PIdata ) ) ;

  // I think this is for the vector  
  if( read_prop( prop , S1 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // read the one above it apart from the last one
  if( read_prop( prop , S1UP ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // copy for the final timeslice
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    equate_spinor_minus( &S1END[x] , &S1[x] ) ;
  }

  // loop the timeslices
  for( t = 0 ; t < ( LT-1 ) ; t++ ) {

    // multiple time source support
    const size_t tshifted = ( t + LT - prop.origin[ ND-1 ] ) % LT ;
    
    // parallel loop with an error flag
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT-2 ) ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      #pragma omp for private(x) schedule(dynamic)
      for( x = 0 ; x < LCU ; x++ ) {
	// do the conserved-local contractions
	contract_conserved_local_site( DATA_AA , DATA_VV , 
				       lat , S1 , S1UP , S1 , S1UP ,
				       GAMMAS , AGMAP , VGMAP , x , 
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
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
      memcpy( &S1UP[x] , &S1f[x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] cl-flavour diagonal done %.f %%", 
	   (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }

  // and contract the final timeslice
  const size_t tshifted = ( t + LT - prop.origin[ ND-1 ] ) % LT ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    contract_conserved_local_site( DATA_AA , DATA_VV , 
				   lat , S1 , S1END , S1 , S1END ,
				   GAMMAS , AGMAP , VGMAP , x , 
				   tshifted ) ;
  }
  printf("\r[VPF] cl-flavour diagonal done 100%% \n" ) ; 

 memfree :

  // free our spinors
  free( S1f ) ; free( S1 ) ; free( S1UP ) ; free( S1END ) ;

  if( error_code != FAILURE ) {

    // derivatives delta_\mu V_\mu(x)
    WI_configspace_bwd( DATA_VV , lat ) ;

    // compute the time moments
    tmoments( DATA_AA , DATA_VV , outfile , CONSERVED_LOCAL ) ;

    // do all the momspace stuff away from the contractions
    if( prop.source == POINT ) { 
      momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		       CONSERVED_LOCAL ) ;
    }
  }
  print_time( ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  return error_code ;
}
