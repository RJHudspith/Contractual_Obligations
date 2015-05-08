/**
   @file bar_contractions.c
   @brief baryon contraction codes
 */

#include "common.h"

#include "bar_ops.h"      // baryon operations
#include "contractions.h" // gamma_mul_lr()

// baryon contraction of ( \Gamma ) ( \Gamma^{T} ) S3 ( S2 S1 )
void
baryon_contract_site( double complex **term ,
		      const struct spinor S1 , 
		      const struct spinor S2 , 
		      const struct spinor S3 , 
		      const struct gamma Cgmu ,
		      const struct gamma CgmuT )
{
  // get diquark
  struct spinor DiQ = S1 ;
  gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;

  // Cross color product and sink Dirac trace back into DiQ
  cross_color_trace( &DiQ , S2 ) ;

  // loop over open dirac indices
  int odc , OD1 , OD2 ;
  for( odc = 0 ; odc < NSNS ; odc++ ) {
    // open dirac index for source and sink
    OD1 = odc / NS ;
    OD2 = odc % NS ;
    // Contract with the final propagator and trace out the source Dirac indices
    // A polarization must still be picked for the two open Dirac indices offline
#if NS == 4
    term[0][ odc ] += baryon_contract( DiQ, S3 , 0 , 0 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 1 , 1 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 2 , 2 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 3 , 3 , OD1 , OD2 ) ;

    term[1][ odc ] += baryon_contract( DiQ, S3 , 0 , OD1 , 0 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 1 , OD1 , 1 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 2 , OD1 , 2 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 3 , OD1 , 3 , OD2 ) ; 
#else
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      term[0][ odc ] += baryon_contract( DiQ, S3 , dirac , dirac , OD1 , OD2 ) ;
      term[1][ odc ] += baryon_contract( DiQ, S3 , dirac , OD1 , dirac , OD2 ) ;
    }
#endif
  }
  return ;
}

// baryon contraction of ( \Gamma ) ( \Gamma^{T} ) S3 ( S2 S1 )
// accumulates site-wise value in flattened "in" array
void
baryon_contract_site_mom( double complex **in ,
			  const struct spinor S1 , 
			  const struct spinor S2 , 
			  const struct spinor S3 , 
			  const struct gamma Cgmu ,
			  const struct gamma CgmuT ,
			  const int GSRC ,
			  const int site )
{
  // get diquark
  struct spinor DiQ = S1 ;
  gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;

  // Cross color product and sink Dirac trace back into DiQ
  cross_color_trace( &DiQ , S2 ) ;

  // loop over open dirac indices
  int odc , OD1 , OD2 ;
  for( odc = 0 ; odc < NSNS ; odc++ ) {
    // open dirac index for source and sink
    OD1 = odc / NS ;
    OD2 = odc % NS ;
    // Contract with the final propagator and trace out the source Dirac indices
    // A polarization must still be picked for the two open Dirac indices offline
    #if NS == 4
    // first term
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 0 , 0 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 1 , 1 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 2 , 2 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 3 , 3 , OD1 , OD2 ) ;
    // second term
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 0 , OD1 , 0 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 1 , OD1 , 1 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 2 , OD1 , 2 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 3 , OD1 , 3 , OD2 ) ;
    #else
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , dirac , dirac , OD1 , OD2 ) ;
      in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , dirac , OD1 , dirac , OD2 ) ;
    }
    #endif
  }
  return ;
}

// baryon contraction for the omega
void
baryon_contract_omega_site( double complex **term ,
			    const struct spinor S1 , 
			    const struct spinor S2 , 
			    const struct spinor S3 , 
			    const struct gamma Cgmu ,
			    const struct gamma CgmuT )
{
  // get diquark
  struct spinor DiQ = S1 , CgS51 = S1 , CgS52 = S1 ;
  gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;
  gamma_mul_l( &CgS51 , CgmuT ) ;
  gamma_mul_r( &CgS52 , Cgmu ) ;

  // Cross color product and sink Dirac trace back into DiQ
  struct spinor DiQ51 = CgS51 , DiQ52 = CgS52 ;
  cross_color_trace( &DiQ , S2 ) ;
  cross_color_trace( &DiQ51 , CgS52 ) ;
  cross_color_trace( &DiQ52 , S2 ) ;

  // loop over open dirac indices
  int odc , OD1 , OD2 ;
  for( odc = 0 ; odc < NSNS ; odc++ ) {
    // open dirac index for source and sink
    OD1 = odc / NS ;
    OD2 = odc % NS ;
    // Contract with the final propagator and trace out the source Dirac indices
    // A polarization must still be picked for the two open Dirac indices offline
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      term[0][odc] += baryon_contract( DiQ , S3 , dirac , dirac , OD1 , OD2 ) ;
      term[1][odc] += baryon_contract( DiQ , S3 , dirac , OD1 , dirac , OD2 ) ;
      term[2][odc] += baryon_contract( DiQ51 , S3 , dirac , dirac , OD1 , OD2 ) ;
      term[3][odc] += baryon_contract( DiQ52 , CgS52, dirac , OD1 , dirac , OD2 ) ;
      term[4][odc] += baryon_contract( DiQ52 , CgS52, OD1 , dirac , dirac , OD2 ) ;
      term[5][odc] += baryon_contract( DiQ51 , S3 , OD1 , dirac , dirac , OD2 ) ;
    }
  }
  return ;
}

