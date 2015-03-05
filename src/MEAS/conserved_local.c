/**
   @file conserved_local.c
   @brief conserved-local Wilson currents

   For the vector ...
  \[ Tr[ ( 1 - \gamma_\mu )U_\mu S( x + \mu ) gamma_\nu adj( S( x ) ) 
      + ( 1 + \gamma_\mu) U_\mu^{\dagger} S( x ) \gamma_\nu \adj( S( x + \mu ) ) ]\]
  
  Doing the four contractions
  \[
    sum += Tr( Id * U_\mu(x) S( x + \mu ) \gamma_\nu adj( S( x ) ) ],\\
    sum += -Tr( gamma_\mu U_\mu(x) S( x + \mu ) \gamma_\nu adj( S( x ) ) ],\\
    sum += Tr( Id * U_\mu^{\dagger}(x) S( x ) \gamma_\nu adj( S( x + \mu ) ) ],\\
    sum += Tr( gamma_\mu U_\mu^{\dagger}(x) S( x ) \gamma_\nu adj( S( x + \mu ) ) ].
   \]
   Once we have U_\mu S( x + \mu ) computed we can ship this out to our meson contracting code right?
 */

#include "common.h"

#include "gammas.h"          // gamma matrices
#include "io.h"              // read_prop
#include "momspace_PImunu.h" // momentum space VPF

// compute the conserved local for a correlator
int
conserved_local( FILE *fprop1 , 
		 const proptype proptype1 ,
		 const struct site *lat )
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
  struct spinor *S1UP = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1END = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // I think this is for the vector  
  read_prop( fprop1 , S1 , proptype1 ) ;

  // copy for the final timeslice
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    int d1 , d2 ;
    // crossing a boundary flips the sign only for antiperiodic BC
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	int c1, c2 ;
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    S1END[x].D[d1][d2].C[c1][c2] = -S1[x].D[d1][d2].C[c1][c2] ;
	  }
	}
	//
      }
    }
  }

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0-1 ; t++ ) {

    // read the one above it apart from the last one
    read_prop( fprop1 , S1UP , proptype1 ) ;

    // do the conserved-local contractions
    contract_conserved_local( DATA_AA , DATA_VV , 
			      lat , S1 , S1UP , S1 , S1UP ,
			      GAMMAS , AGMAP , VGMAP , t ) ;

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
    }
  }

  // and contract the final timeslice
  contract_conserved_local( DATA_AA , DATA_VV , 
			    lat , S1 , S1END , S1 , S1END ,
			    GAMMAS , AGMAP , VGMAP , t ) ;

  // free our spinors
  free( S1 ) ;
  free( S1UP ) ;
  free( S1END ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // do all the momspace stuff away from the contractions
  momspace_PImunu( DATA_AA , DATA_VV ) ;

  // free the AA & VV data
  free( DATA_AA ) ;
  free( DATA_VV ) ;

  return SUCCESS ;
}
