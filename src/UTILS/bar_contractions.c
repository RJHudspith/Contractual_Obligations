/**
   @file bar_contractions.c
   @brief baryon contraction codes
 */

#include "common.h"

#include "bar_ops.h"      // baryon operations
#include "contractions.h" // gamma_mul_lr()
#include "gammas.h"       // Cgmu and CgmuT

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

    term[1][ odc ] += baryon_contract( DiQ, S3 , OD1 , 0 , 0 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , OD1 , 1 , 1 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , OD1 , 2 , 2 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , OD1 , 3 , 3 , OD2 ) ; 
#else
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      term[0][ odc ] += baryon_contract( DiQ, S3 , dirac , dirac , OD1 , OD2 ) ;
      term[1][ odc ] += baryon_contract( DiQ, S3 , OD1 , dirac , dirac , OD2 ) ;
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
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 0 , 0 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 1 , 1 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 2 , 2 , OD1 , OD2 ) ;
    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , 3 , 3 , OD1 , OD2 ) ;

    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , OD1 , 0 , 0 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , OD1 , 1 , 1 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , OD1 , 2 , 2 , OD2 ) ;
    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , OD1 , 3 , 3 , OD2 ) ;

    #else
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , dirac , dirac , OD1 , OD2 ) ;
      in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] += baryon_contract( DiQ, S3 , OD1 , dirac , dirac , OD2 ) ;
    }
    #endif
  }
  return ;
}

// must be called within a parallel environment
void
baryon_contract_walls( struct mcorr **Buud_corrWW , 
		       struct mcorr **Buuu_corrWW ,
		       struct mcorr **Buds_corrWW ,
		       const struct spinor SUM1 ,
		       const struct spinor SUM2 ,
		       const struct spinor SUM3 ,
		       const struct gamma *GAMMAS ,
		       const int t )
{
  // accumulate the sums with open dirac indices
  int GSRC ;
#pragma omp for private(GSRC) schedule(dynamic)
  for( GSRC = 0 ; GSRC < ( B_CHANNELS ) ; GSRC++ ) {
    // allocate these locally
    double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
    term[0] = calloc( NSNS , sizeof( double complex ) ) ;
    term[1] = calloc( NSNS , sizeof( double complex ) ) ;
    // recompute teh Cgmus, these are basically free to calculate
    const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
    const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;
    baryon_contract_site( term , SUM1 , SUM1 , SUM1 , Cgmu , CgmuT ) ;
    // wall contractions project to zero spatial momentum explicitly
    int odc ;
    for( odc = 0 ; odc < NSNS ; odc++ ) {
      Buds_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = term[0][odc] ;
      Buud_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = term[0][odc] + term[1][odc] ;
      Buuu_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = 2 * term[0][odc] + 4 * term[1][odc] ;
    }
    free( term[0] ) ; free( term[1] ) ; free( term ) ;
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
      term[1][odc] += baryon_contract( DiQ , S3 , OD1 , dirac , dirac , OD2 ) ;
      term[2][odc] += baryon_contract( DiQ51 , S3 , dirac , dirac , OD1 , OD2 ) ;
      term[3][odc] += baryon_contract( DiQ52 , CgS52, OD1 , dirac , dirac , OD2 ) ;
      term[4][odc] += baryon_contract( DiQ52 , CgS52, dirac , OD1 , dirac , OD2 ) ;
      term[5][odc] += baryon_contract( DiQ51 , S3 , dirac , OD1 , dirac , OD2 ) ;
    }
  }
  return ;
}

// baryonic momentum project
void
baryon_momentum_project( struct mcorr **Buud_corr , 
			 struct mcorr **Buuu_corr ,
			 struct mcorr **Buds_corr ,
			 double complex **in ,
			 double complex **out ,
			 const void *forward ,
			 const void *backward ,
			 const struct veclist *list ,
			 const int NMOM[ 1 ] ,
			 const int t )
{
  int GSodc ; // flatteded open dirac indices
#ifdef HAVE_FFTW3_H
  const fftw_plan *forw = ( const fftw_plan* )forward ;
  // perform the FFTS separately here
  #pragma omp parallel for private(GSodc)
  for( GSodc = 0 ; GSodc < ( B_CHANNELS * NSNS ) ; GSodc++ ) {
    const int GSRC = GSodc / NSNS ;
    const int odc = GSodc % NSNS ;
    const int idx = 2 * GSodc ;
    fftw_execute( forw[ 0 + idx ] ) ;
    fftw_execute( forw[ 1 + idx ] ) ;
    const double complex *sum1 = out[ 0 + idx ] ;
    const double complex *sum2 = out[ 1 + idx ] ;
    int p ;
    for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
      const int lid = list[ p ].idx ;
      Buds_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = sum1[ lid ] ;
      Buud_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = sum1[ lid ] + sum2[ lid ] ;
      Buuu_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = 2 * sum1[ lid ] + 4 * sum2[ lid ] ;
    }
  }
#else
  #pragma omp parallel for private(GSodc)
  for( GSodc = 0 ; GSodc < ( B_CHANNELS * NSNS ) ; GSodc++ ) {
    const int GSRC = GSodc / NSNS ;
    const int odc = GSodc % NSNS ;
    const int idx = 2 * GSodc ;
    register double complex sum1 = 0.0 , sum2 = 0.0 ;
    int site ;
    for( site = 0 ; site < LCU ; site++ ) {
      sum1 += in[ 0 + idx ][ site ] ;
      sum2 += in[ 1 + idx ][ site ] ;
    }
    Buds_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = sum1 ;
    Buud_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = sum1 + sum2 ;
    Buuu_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = 2 * sum1 + 4 * sum2 ;
  }
#endif
  return ;
}
