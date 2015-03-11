/**
   @file local_local.c
   @brief local_local current contractions
 */

#include "common.h"

#include "gammas.h"          // gamma matrices
#include "io.h"              // read_prop
#include "matrix_ops.h"      // constant_mul_gauge
#include "momspace_PImunu.h" // momentum space VPF

// compute the conserved local for a correlator
int
local_local( FILE *prop1 , 
	     const proptype proptype1 ,
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
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {

    // read the one above it apart from the last one
    if( read_prop( prop1 , S1 , proptype1 ) == FAILURE ) {
      free( GAMMAS ) ;
      free( S1 ) ;
      free( DATA_AA ) ;
      free( DATA_VV ) ;
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
  momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		   LOCAL_LOCAL ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop1 ) ;
  read_check_header( prop1 , GLU_FALSE ) ;

  return SUCCESS ;
}

// do the local-local for different fermion props
int
local_local_double( FILE *prop1 , 
		    const proptype proptype1 ,
		    FILE *prop2 , 
		    const proptype proptype2 ,
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
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S2 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {

    // read the propagators for this timeslice
    if( read_prop( prop1 , S1 , proptype1 ) || 
	read_prop( prop1 , S2 , proptype2 ) == FAILURE ) {
      free( GAMMAS ) ;
      free( S1 ) ;
      free( DATA_AA ) ;
      free( DATA_VV ) ;
      return FAILURE ;
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
  momspace_PImunu( DATA_AA , DATA_VV , CUTINFO , outfile ,
		   LOCAL_LOCAL ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  // rewind file and read header again
  rewind( prop1 ) ;
  read_check_header( prop1 , GLU_FALSE ) ;
  rewind( prop2 ) ;
  read_check_header( prop2 , GLU_FALSE ) ;

  return SUCCESS ;
}
