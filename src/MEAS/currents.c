/**
   @file currents.c
   @brief various current implementations and contractions
 */

#include "common.h"

#include "contractions.h" // meson_contract
#include "spinor_ops.h"   // spinor - color matrix multiply

/**
   @brief Wilson non-local non-conserved axial
   Computes \f[ 
   \text{Tr}\left[ \gamma_\nu \gamma_5 ( \gamma_5  S2^{\dagger}(x) \gamma_5 ) \gamma_\mu \gamma_5 U_\mu\left( x+ a\frac{\hat\mu}{2} \right)S1(x+\mu) \right] + \text{Tr}\left[ \gamma_\nu \gamma_5 ( \gamma_5  S2^{\dagger}(x+\mu) \gamma_5 ) \gamma_\mu \gamma_5 U_\mu\left( x+ a\frac{\hat\mu}{2} \right)^{\dagger} S1(x) \right] 
   \f]
*/
static double complex
CL_munu_AA( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	    const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	    const struct spinor S2 ,       // S_2
	    const struct spinor S2xpmu ,   // S_2( x + \mu )
	    const struct gamma *GAMMAS ,
	    const int mu ,
	    const int nu )
{
  return \
    0.5 * ( meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ mu ] , US1xpmu , GAMMAS[ GAMMA_5 ] ) +
	    meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ mu ] , UdS1x , GAMMAS[ GAMMA_5 ] ) ) ;
}

/**
   @brief Wilson non-local non-conserved axial
   Computes \f[ 
   \text{Tr}\left[ \gamma_\nu ( \gamma_5  S2^{\dagger}(x) \gamma_5 ) \gamma_\mu U_\mu\left( x+ a\frac{\hat\mu}{2} \right)S1(x+\mu) \right] + \text{Tr}\left[ \gamma_\nu ( \gamma_5  S2^{\dagger}(x+\mu) \gamma_5 ) \gamma_\mu U_\mu\left( x+ a\frac{\hat\mu}{2} \right)^{\dagger} S1(x) \right]    \f]
*/
// non-conserved non-local vector current 
static double complex
NCL_munu_VV( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	     const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	     const struct spinor S2 ,       // S_2
	     const struct spinor S2xpmu ,   // S_2( x + \mu )
	     const struct gamma *GAMMAS ,
	     const int mu ,
 	     const int nu )
{
  return \
    0.5 * ( meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ mu ] , US1xpmu , GAMMAS[ GAMMA_5 ] ) +
	    meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ mu ] , UdS1x , GAMMAS[ GAMMA_5 ] ) ) ;
}

// Conserved-Local Vector current
static double complex
CL_munu_VV( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	    const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	    const struct spinor S2 ,       // S_2
	    const struct spinor S2xpmu ,   // S_2( x + \mu )
	    const struct gamma *GAMMAS ,
	    const int mu ,
	    const int nu )
{
  return \
    0.5 * (
	   -meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ IDENTITY ] , US1xpmu , GAMMAS[ GAMMA_5 ] )
	   +meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ mu ] , US1xpmu , GAMMAS[ GAMMA_5 ] )
	   +meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ IDENTITY ], UdS1x , GAMMAS[ GAMMA_5 ] )
	   +meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ mu ] , UdS1x , GAMMAS[ GAMMA_5 ] )
	    ) ;
}

// man this has a lot of arguments -> TODO :: reduce these somehow
void
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
  int x ;
#pragma omp parallel for private(x)
  for( x = 0 ; x < LCU ; x++ ) {

    struct gamma G1 , G2 ;

    struct spinor US1xpmu , UdS1x , S2xpmu ; // temporary storage for the gauge-multiplied
    const int i = x + LCU * t ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {

      // multiply through by the link matrix 
      // if we are in the t-direction we use the "UP" space
      if( mu == ND-1 ) {
	gauge_spinor( &US1xpmu , lat[i].O[mu] , S1UP[ x ] ) ;  // U S(x+\mu)
	S2xpmu = S2UP[ x ] ;
      } else {
	const int xpmu = lat[x].neighbor[mu] ;
	gauge_spinor( &US1xpmu , lat[i].O[mu] , S1[ xpmu ] ) ; // U S(x+\mu)
	S2xpmu = S2[ xpmu ] ;
      }
      gaugedag_spinor( &UdS1x , lat[i].O[mu] , S1[x] ) ;       // U^{\dagger} S(x)

      for( nu = 0 ; nu < ND ; nu++ ) {

	// I need to think about the axial
	DATA_AA[i].PI[mu][nu] = CL_munu_AA( US1xpmu , UdS1x , 
					    S1[x] , S2xpmu ,
					    GAMMAS , 
					    AGMAP[ mu ] , AGMAP[ nu ] ) ;
	
	// vectors 
	DATA_VV[i].PI[mu][nu] = CL_munu_VV( US1xpmu , UdS1x , 
					    S1[x] , S2xpmu ,
					    GAMMAS ,
					    VGMAP[mu] , VGMAP[ nu ] ) ;
      }
    }
  }
  return ;
}

// LL
void
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

	DATA_AA[i].PI[mu][nu] = \
	  meson_contract( GAMMAS[ AGMAP[ nu ] ] , S2[x] , 
			  GAMMAS[ AGMAP[ mu ] ] , S1[x] ,
			  GAMMAS[ GAMMA_5 ] ) ;

	DATA_VV[i].PI[mu][nu] = \
	  meson_contract( GAMMAS[ VGMAP[ nu ] ] , S2[x] , 
			  GAMMAS[ VGMAP[ mu ] ] , S1[x] ,
			  GAMMAS[ GAMMA_5 ] ) ;
      }
    }
  }
  return ;
}
