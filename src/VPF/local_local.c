/**
   @file local_local.c
   @brief local_local current contractions
 */

#include "common.h"

#include "basis_conversions.h" // rotate slice
#include "currents.h"          // contractions
#include "gammas.h"            // gamma matrices
#include "io.h"                // read_prop
#include "matrix_ops.h"        // constant_mul_gauge
#include "momspace_PImunu.h"   // momentum space VPF
#include "tmoments_PImunu.h"   // time moments calculations
#include "read_propheader.h"   // read_propheader()

// compute the conserved local for a correlator
int
local_local( struct propagator prop ,
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
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // and our spinor
  struct spinor *S1 ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free( S1 ) ; free( GAMMAS ) ;
    printf( "[MESONS] memalign failure \n" ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {

    // read the one above it apart from the last one
    if( read_prop( prop , S1 ) == FAILURE ) {
      free( GAMMAS ) ; free( S1 ) ;
      free( DATA_AA ) ; free( DATA_VV ) ;
      return FAILURE ;
    }

    // do the conserved-local contractions
    contract_local_local( DATA_AA , DATA_VV , S1 , S1 , 
			  GAMMAS , AGMAP , VGMAP , t ) ;

    // status
    printf("\r[VPF] ll-flavour diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // free our spinor(s)
  free( S1 ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // temporal data too
  tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  if( prop.source == POINT ) {
    momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		     LOCAL_LOCAL ) ;
  }

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  return SUCCESS ;
}

// do the local-local for different fermion props
int
local_local_double( struct propagator prop1 ,
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

  // and our spinors
  struct spinor *S1 ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free( S1 ) ; free( GAMMAS ) ;
    printf( "[MESONS] memalign failure \n" ) ;
    return FAILURE ;
  }
  struct spinor *S2 ;
  if( posix_memalign( (void**)&S2 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free( S1 ) ; free( S2 ) ; free( GAMMAS ) ;
    printf( "[MESONS] memalign failure \n" ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {

    // read the propagators for this timeslice
    if( read_prop( prop1 , S1 ) || 
	read_prop( prop2 , S2 ) == FAILURE ) {
      free( GAMMAS ) ; free( S1 ) ;
      free( DATA_AA ) ; free( DATA_VV ) ;
      return FAILURE ;
    }

    // rotate if necessary?
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
      nrel_rotate_slice( S2 ) ;
    }

    // do the conserved-local contractions
    contract_local_local( DATA_AA , DATA_VV , S1 , S2 , 
			  GAMMAS , AGMAP , VGMAP , t ) ;

    // status
    printf("\r[VPF] ll-flavour off diagonal done %.f %%", 
	   (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // free our spinor(s)
  free( S1 ) ;
  free( S2 ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // time moment data
  tmoments( DATA_AA , DATA_VV , outfile , LOCAL_LOCAL ) ;

  // do all the momspace stuff away from the contractions
  if( prop1.source == POINT ) { 
    momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		     LOCAL_LOCAL ) ;
  }

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  return SUCCESS ;
}
