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

#include "correlators.h"     // our contraction code
#include "gammas.h"          // gamma matrices
#include "io.h"              // read_prop
#include "matrix_ops.h"      // color matrix multiply
#include "momspace_PImunu.h" // momentum space VPF

// grim
static double complex
contract_munu( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	       const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	       const struct spinor S2 ,       // S_2
	       const struct spinor S2xpmu ,   // S_2( x + \mu )
	       const struct gamma *GAMMAS ,
	       const int mu ,
	       const int nu )
{
  register double complex sum = 0.0 ;
  sum += local_meson_correlator( US1xpmu , S2 , GAMMAS[ GAMMA_5 ] , GAMMAS[ IDENTITY ] , GAMMAS[ nu ] ) ;
  sum -= local_meson_correlator( US1xpmu , S2 , GAMMAS[ GAMMA_5 ] , GAMMAS[ mu ] , GAMMAS[ nu ] ) ;
  sum += local_meson_correlator( UdS1x , S2xpmu , GAMMAS[ GAMMA_5 ] , GAMMAS[ IDENTITY ] , GAMMAS[ nu ] ) ;
  sum += local_meson_correlator( UdS1x , S2xpmu , GAMMAS[ GAMMA_5 ] , GAMMAS[ mu ] , GAMMAS[ nu ] ) ;
  return sum * 0.5 ;
}
	       
// man this has a lot of arguments -> TODO :: reduce these somehow
static void
contract_conserved_local( struct PIdata *DATA_AA ,
			  struct PIdata *DATA_VV ,
			  const struct site *lat ,
			  const struct spinor *S1 ,
 			  const struct spinor *S1UP ,
			  const struct spinor *S2 ,
			  const struct spinor *S2UP ,
			  const struct gamma *GAMMAS ,
			  const int AGMAP[ ND ] ,
			  const int VGMAP[ ND ] ,
			  const int t ) 
{
  // List of obvious optimisations so I don't forget
  // 
  // the correlator has to be real -> halves the time of the contraction
  // axial is basically the same after VGMAP -> AGMAP
  int x ;
  #pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    struct spinor US1xpmu , UdS1x , S2xpmu ; // temporary storage for the gauge-multiplied
    const int i = x + LCU * t ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {

      // multiply through by the link matrix 

      // if we are in the t-direction we use the "UP" space
      if( mu == ND-1 ) {
	gauge_spinor( &US1xpmu , S1UP[ x ] , lat[i].O[mu] ) ;  // U S(x+\mu)
	gaugedag_spinor( &UdS1x , S1[ x ] , lat[i].O[mu] ) ;   // U^{\dagger} S(x)
	S2xpmu = S1UP[ x ] ;
      } else {
	const int xpmu = lat[x].neighbor[mu] ;
	gauge_spinor( &US1xpmu , S1[ xpmu ] , lat[i].O[mu] ) ; // U S(x+\mu)
	gaugedag_spinor( &UdS1x , S1[ x ] , lat[i].O[mu] ) ;   // U^{\dagger} S(x)
	S2xpmu = S1[ xpmu ] ;
      }

      for( nu = 0 ; nu < ND ; nu++ ) {
	// I need to think about the axial

	// vectors 
	DATA_VV[i].PI[mu][nu] = contract_munu( US1xpmu , UdS1x , 
					       S1[x] , S2xpmu ,
					       GAMMAS , VGMAP[ mu ] , VGMAP[ nu ] ) ;
      }
    }
  }

  return ;
}

// LL
static void
contract_local_local( struct PIdata *DATA_AA ,
		      struct PIdata *DATA_VV ,
		      const struct spinor *S1 ,
		      const struct spinor *S2 ,
		      const struct gamma *GAMMAS ,
		      const int AGMAP[ ND ] ,
		      const int VGMAP[ ND ] ,
		      const int t ) 
{
  int x ;
  #pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    const int i = x + LCU * t ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	// axial-axial
	DATA_AA[i].PI[mu][nu] = local_meson_correlator( S1[x] , S2[x] , GAMMAS[ GAMMA_5 ] , 
							GAMMAS[ AGMAP[ mu ] ] , GAMMAS[ AGMAP[ nu ] ] ) ;

	// vector-vector
	DATA_VV[i].PI[mu][nu] = local_meson_correlator( S1[x] , S2[x] , GAMMAS[ GAMMA_5 ] , 
							GAMMAS[ VGMAP[ mu ] ] , GAMMAS[ VGMAP[ nu ] ] ) ;
      }
    }
  }
  return ;
}

// compute the conserved local for a correlator
int
conserved_local( FILE *fprop1 , 
		 const proptype proptype1 ,
		 const struct site *lat )
{
  const int VGMAP[ ND ] = { GAMMA_0 , GAMMA_1 , GAMMA_2 , GAMMA_3 } ;
#ifdef CHROMA_DIRAC_CONVENTION
  const int AGMAP[ ND ] = { GAMMA_5^GAMMA_0 , GAMMA_5^GAMMA_1 , GAMMA_5^GAMMA_2 , GAMMA_5^GAMMA_3 } ;
  //printf( "CHROMA convention unsupported ... Leaving \n" ) ;
  //return FAILURE ;
#else
  // need to look these up
  const int AGMAP[ ND ] = { GAMMA_5 + 1 , GAMMA_5 + 2 , GAMMA_5 + 3 , GAMMA_5 + 4 } ;
#endif

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS , proptype1 ) ;

  // over the whole volume is not as expensive as you may think
  struct PIdata *DATA_AA = malloc( LVOLUME * sizeof( struct PIdata ) ) ;
  struct PIdata *DATA_VV = malloc( LVOLUME * sizeof( struct PIdata ) ) ;

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S1UP = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S1END = malloc( VOL3 * sizeof( struct spinor ) ) ;

  // I think this is for the vector  
  read_prop( fprop1 , S1 , proptype1 ) ;

  // copy for the final timeslice
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {
    memcpy( &S1END[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
  }

  // loop the timeslices
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {

    // read the one above it apart from the last one
    if( t != L0-1 ) {
      read_prop( fprop1 , S1UP , proptype1 ) ;
    } else {
      //
      #pragma omp parallel for private(x)
      for( x = 0 ; x < LCU ; x++ ) {
	memcpy( &S1UP[x] , &S1END[x] , sizeof( struct spinor ) ) ;
      }
      //
    }

    // do the contractions
    //#if 0
    contract_conserved_local( DATA_AA , DATA_VV , 
			      lat , S1 , S1UP , S1 , S1UP ,
			      GAMMAS , AGMAP , VGMAP , t ) ;
    // #endif
    #if 0
    contract_local_local( DATA_AA , DATA_VV ,
	  		  S1 , S1 ,
			  GAMMAS , AGMAP , VGMAP , t ) ;
    #endif

    #pragma omp parallel for private(x)
    for( x = 0 ; x < LCU ; x++ ) {
      // copy spinors over a timeslice
      memcpy( &S1[x] , &S1UP[x] , sizeof( struct spinor ) ) ;
    }
  }

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
