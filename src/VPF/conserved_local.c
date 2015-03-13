/**
   @file conserved_local.c
   @brief conserved-local Wilson currents

   At the moment I treat the NREL like the chiral which is
   probably wrong, but we cross that bridge when we get to it
 */

#include "common.h"

#include "gammas.h"          // gamma matrices
#include "io.h"              // read_prop
#include "matrix_ops.h"      // constant_mul_gauge
#include "momspace_PImunu.h" // momentum space VPF
#include "tmoments_PImunu.h" // time moments too

// compute the conserved local for a correlator
int
conserved_local( struct propagator prop ,
		 const struct site *lat ,
		 const struct cut_info CUTINFO ,
		 const char *outfile )
{
  // vector gamma map
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

#ifdef CHROMA_DIRAC_CONVENTION
  // two of these have the wrong sign this needs to be fixed!!
  const int AGMAP[ ND ] = { GAMMA_5^GAMMA_0 , GAMMA_5^GAMMA_1 , GAMMA_5^GAMMA_2 , GAMMA_5^GAMMA_3 } ;
#else
  // need to look these up
  const int AGMAP[ ND ] = { 6 , 7 , 8 , 9 } ;
#endif

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1UP = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1END = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // I think this is for the vector  
  if( read_prop( prop , S1 ) == FAILURE ) {
    free( S1 ) ; free( S1UP ) ; free( S1END ) ;
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // copy for the final timeslice
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    int d1 , d2 ;
    // crossing a boundary flips the sign only for antiperiodic BC
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	constant_mul_gauge( (double complex*)S1END[x].D[d1][d2].C , 
			    -1 , (double complex*)S1[x].D[d1][d2].C ) ;
      }
    }
  }

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0-1 ; t++ ) {

    // read the one above it apart from the last one
    if( read_prop( prop , S1UP ) == FAILURE ) {
      free( S1 ) ; free( S1UP ) ; free( S1END ) ;
      free( DATA_AA ) ; free( DATA_VV ) ;
      free( GAMMAS ) ;
      return FAILURE ;
    }

    // do the conserved-local contractions
    contract_conserved_local( DATA_AA , DATA_VV , 
			      lat , S1 , S1UP , S1 , S1UP ,
			      GAMMAS , AGMAP , VGMAP , t ) ;

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] cl-flavour diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }

  // and contract the final timeslice
  contract_conserved_local( DATA_AA , DATA_VV , 
			    lat , S1 , S1END , S1 , S1END ,
			    GAMMAS , AGMAP , VGMAP , t ) ;

  printf("\r[VPF] cl-flavour diagonal done 100%% \n" ) ; 

  // derivatives delta_\mu V_\mu(x)
  WI_configspace( DATA_VV , lat ) ;

  // free our spinors
  free( S1 ) ;
  free( S1UP ) ;
  free( S1END ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  tmoments( DATA_AA , DATA_VV , outfile , CONSERVED_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		   CONSERVED_LOCAL ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_check_header( &prop , GLU_FALSE ) ;

  return SUCCESS ;
}

// compute the conserved local for a correlator
int
conserved_local_double( struct propagator prop1 ,
			struct propagator prop2 ,
			const struct site *lat ,
			const struct cut_info CUTINFO ,
			const char *outfile )
{
  // vector gamma map
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;

#ifdef CHROMA_DIRAC_CONVENTION
  // two of these have the wrong sign this needs to be fixed!!
  const int AGMAP[ ND ] = { GAMMA_5^GAMMA_0 , GAMMA_5^GAMMA_1 , GAMMA_5^GAMMA_2 , GAMMA_5^GAMMA_3 } ;
#else
  // need to look these up
  const int AGMAP[ ND ] = { GAMMA_5 + 1 , GAMMA_5 + 2 , GAMMA_5 + 3 , GAMMA_5 + 4 } ;
#endif

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // and our spinor(s)
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1UP = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1END = calloc( VOL3 , sizeof( struct spinor ) ) ;

  struct spinor *S2 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S2UP = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S2END = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // free everything if the read fails
  if( read_prop( prop1 , S1 ) == FAILURE || 
      read_prop( prop2 , S2 ) == FAILURE ) {
    free( S1 ) ; free( S1UP ) ; free( S1END ) ;
    free( S2 ) ; free( S2UP ) ; free( S2END ) ;
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // if we are doing nonrel-chiral mesons we switch chiral to nrel
  if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
    nrel_rotate_slice( S1 ) ;
  } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
    nrel_rotate_slice( S2 ) ;
  }

  // copy for the final timeslice
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    int d1 , d2 ;
    // crossing a boundary flips the sign only for antiperiodic BC
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	constant_mul_gauge( (double complex*)S1END[x].D[d1][d2].C , 
			    -1 , (double complex*)S1[x].D[d1][d2].C ) ;
	constant_mul_gauge( (double complex*)S2END[x].D[d1][d2].C , 
			    -1 , (double complex*)S2[x].D[d1][d2].C ) ;
      }
    }
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0-1 ; t++ ) {

    // read the one above it apart from the last one
    if( read_prop( prop1 , S1UP ) == FAILURE ||
	read_prop( prop2 , S2UP ) == FAILURE ) {
      free( S1 ) ; free( S1UP ) ; free( S1END ) ;
      free( S2 ) ; free( S2UP ) ; free( S2END ) ;
      free( DATA_AA ) ; free( DATA_VV ) ;
      free( GAMMAS ) ;
      return FAILURE ;
    }

    // convert to NREL if needed?
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1UP ) ;
    } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
      nrel_rotate_slice( S2UP ) ;
    }

    // do the conserved-local contractions
    contract_conserved_local( DATA_AA , DATA_VV , 
			      lat , S1 , S1UP , S2 , S2UP ,
			      GAMMAS , AGMAP , VGMAP , t ) ;

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
      memcpy( &S2[x] , &S2UP[x] , sizeof( struct spinor ) ) ;
    }

    // status
    printf("\r[VPF] cl-flavour diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }

  // and contract the final timeslice
  contract_conserved_local( DATA_AA , DATA_VV , 
			    lat , S1 , S1END , S2 , S2END ,
			    GAMMAS , AGMAP , VGMAP , t ) ;

  printf("\r[VPF] cl-flavour off diagonal done 100%% \n" ) ; 

  // derivatives delta_\mu V_\mu(x)
  WI_configspace( DATA_VV , lat ) ;

  // free our spinors
  free( S1 ) ;
  free( S1UP ) ;
  free( S1END ) ;

  free( S2 ) ;
  free( S2UP ) ;
  free( S2END ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // time moments are interesting also
  tmoments( DATA_AA , DATA_VV , outfile , CONSERVED_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		   CONSERVED_LOCAL ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_check_header( &prop1 , GLU_FALSE ) ;
  rewind( prop2.file ) ; read_check_header( &prop2 , GLU_FALSE ) ;

  return SUCCESS ;
}
