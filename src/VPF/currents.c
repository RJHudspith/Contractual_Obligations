/**
   @file currents.c
   @brief various current implementations and contractions
 */

#include "common.h"

#include "contractions.h" // meson_contract
#include "gammas.h"       // gt_Gdag_gt
#include "matrix_ops.h"   // link multiplies
#include "spinor_ops.h"   // spinor - color matrix multiply

// non-conserved, non-local Axial current
double complex
CL_munu_AA( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	    const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	    const struct spinor S2 ,       // S_2
	    const struct spinor S2xpmu ,   // S_2( x + \mu )
	    const struct gamma *GAMMAS ,
	    const size_t mu ,
	    const size_t nu )
{
  return \
    0.5 * ( meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ mu ] , US1xpmu , GAMMAS[ GAMMA_5 ] ) +
	    meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ mu ] , UdS1x , GAMMAS[ GAMMA_5 ] ) ) ;
}

// non-conserved non-local vector current 
double complex
NCL_munu_VV( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	     const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	     const struct spinor S2 ,       // S_2
	     const struct spinor S2xpmu ,   // S_2( x + \mu )
	     const struct gamma *GAMMAS ,
	     const size_t mu ,
 	     const size_t nu )
{
  return \
    0.5 * ( meson_contract( GAMMAS[ nu ] , S2 , GAMMAS[ mu ] , US1xpmu , GAMMAS[ GAMMA_5 ] ) +
	    meson_contract( GAMMAS[ nu ] , S2xpmu , GAMMAS[ mu ] , UdS1x , GAMMAS[ GAMMA_5 ] ) ) ;
}

// Conserved-Local Vector current
double complex
CL_munu_VV( const struct spinor US1xpmu ,  // U S_1( x + \mu )
	    const struct spinor UdS1x ,    // U^{\dagger} S_1( x )
	    const struct spinor S2xpmu ,   // S_2( x + \mu )
	    const struct spinor S2 ,       // S_2
	    const struct gamma *GAMMAS ,
	    const size_t mu ,
	    const size_t nu )
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
contract_conserved_local_site( struct PIdata *DATA_AA ,
			       struct PIdata *DATA_VV ,
			       const struct site *lat ,
			       const struct spinor *S1 ,
			       const struct spinor *S1UP ,
			       const struct spinor *S2 ,
			       const struct spinor *S2UP ,
			       const struct gamma *GAMMAS ,
			       const size_t AGMAP[ ND ] ,
			       const size_t VGMAP[ ND ] ,
			       const size_t x ,
			       const size_t t ) 
{
  struct spinor US1xpmu , UdS1x , S2xpmu ; // temporary storage for the gauge-multiplied
  const size_t i = x + LCU * t ;
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {

    // multiply through by the link matrix 
    // if we are in the t-direction we use the "UP" space
    if( mu == ND-1 ) {
      gauge_spinor( &US1xpmu , lat[i].O[mu] , S1UP[x] ) ;  // U S(x+\mu)
      S2xpmu = S2UP[ x ] ;
    } else {
      const size_t xpmu = lat[x].neighbor[mu] ;
      gauge_spinor( &US1xpmu , lat[i].O[mu] , S1[xpmu] ) ; // U S(x+\mu)
      S2xpmu = S2[ xpmu ] ;
    }
    gaugedag_spinor( &UdS1x , lat[i].O[mu] , S1[x] ) ; // U^{\dagger} S(x)

    for( nu = 0 ; nu < ND ; nu++ ) {

      // I need to think about the axial
      DATA_AA[i].PI[mu][nu] = CL_munu_AA( US1xpmu , UdS1x , 
					  S2xpmu , S2[ x ] , 
					  GAMMAS , 
					  AGMAP[ mu ] , AGMAP[ nu ] ) ;
	
      // vectors 
      DATA_VV[i].PI[mu][nu] = -CL_munu_VV( US1xpmu , UdS1x , 
					   S2xpmu , S2[ x ] ,
					   GAMMAS ,
					   VGMAP[ mu ] , VGMAP[ nu ] ) ;
      //
    }
  }
  //exit(1) ;

  return ;
}

// LL
void
contract_local_local_site( struct PIdata *DATA_AA ,
			   struct PIdata *DATA_VV ,
			   const struct spinor *S1 ,
			   const struct spinor *S2 ,
			   const struct gamma *GAMMAS ,
			   const size_t AGMAP[ ND ] ,
			   const size_t VGMAP[ ND ] ,
			   const size_t x ,
			   const size_t t ) 
{
  const size_t i = x + LCU * t ;
  size_t munu ;
  for( munu = 0 ; munu < ND*ND ; munu++ ) {
    const size_t mu = munu / ND ;
    const size_t nu = munu % ND ;
    DATA_AA[i].PI[mu][nu] =				\
      meson_contract( GAMMAS[ AGMAP[ nu ] ] , S2[ x ] , 
		      GAMMAS[ AGMAP[ mu ] ] , S1[ x ] ,
		      GAMMAS[ GAMMA_5 ] ) ;
    
    DATA_VV[i].PI[mu][nu] =				\
      meson_contract( GAMMAS[ VGMAP[ nu ] ] , S2[ x ] , 
		      GAMMAS[ VGMAP[ mu ] ] , S1[ x ] ,
		      GAMMAS[ GAMMA_5 ] ) ;
  }
  return ;
}
