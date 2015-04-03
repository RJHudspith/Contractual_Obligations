/**
   @file baryons.c
   @brief baryon contraction code


   This needs some splainin
*/

#include "common.h"

#include "bar_contractions.h"  // cross_color trace and baryon contract
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // gamma_mul_lr()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas() && gamma_mmul
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()

// contracts S2 ( S1 C S1 C ) , odd one out is outside
int
baryons_2diagonal( struct propagator prop1 ,
		   struct propagator prop2 ,
		   const char *outfile )
{
  // Define our output correlators, with 6 channels and 16 components
  struct correlator **Buds_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  struct correlator **Buud_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  struct correlator **Buuu_corr = allocate_corrs( B_CHANNELS , NSNS ) ;

  // and our spinors
  struct spinor *S1 = NULL ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free_corrs( Buds_corr , B_CHANNELS , NSNS ) ; 
    free_corrs( Buud_corr , B_CHANNELS , NSNS ) ; 
    free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ; 
    free( S1 ) ;
    printf( "[BARYONS] memalign failure \n" ) ;
    return FAILURE ;
  }
  struct spinor *S2 = NULL ;
  if( posix_memalign( (void**)&S2 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free_corrs( Buds_corr , B_CHANNELS , NSNS ) ; 
    free_corrs( Buud_corr , B_CHANNELS , NSNS ) ; 
    free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ; 
    free( S1 ) ; free( S2 ) ;
    printf( "[BARYONS] memalign failure \n" ) ;
    return FAILURE ;
  }

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( prop1.basis == NREL || prop2.basis == NREL ) { 
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      free( GAMMAS ) ;
      return FAILURE ;
    }
  } else {
    if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      free( GAMMAS ) ;
      return FAILURE ;
    }
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 ) == FAILURE ) {
      free_corrs( Buds_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buud_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ; 
      free( S1 ) ; free( S2 ) ;
      return FAILURE ;
    }
    if( read_prop( prop2 , S2 ) == FAILURE ) {
      free_corrs( Buds_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buud_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ; 
      free( S1 ) ; free( S2 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } else if( prop2.basis == CHIRAL && prop1.basis == NREL ) {
      nrel_rotate_slice( S2 ) ;
    } 

    int GSRC = 0 ;
    // parallelise the furthest out loop
#pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {

      // Define some intermediate spinors
      struct spinor DiQ ;

      // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
      const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
      // precompute \gamma_t ( Cg_\mu )^{*} \gamma_t -> \Gamma^{T} in note
      const struct gamma CgmuT = CGmuT( GAMMAS[ GSRC ] , GAMMAS ) ;

      // accumulate the sums with open dirac indices
      double complex Buds[ NSNS ] = {} ;
      double complex Buud[ NSNS ] = {} ;
      double complex Buuu[ NSNS ] = {} ;

      // loop spatial hypercube
      int site ;
      for( site = 0 ; site < VOL3 ; site++ ) {

	// multiply the di-quark by CgmuT from the left and Cgmu from the right
	DiQ = S1[ site ] ;
	gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;

	// Cross color product and sink Dirac trace back into DiQ
	cross_color_trace( &DiQ , S1[ site ] ) ;

	// loop over open dirac indices
	int odc ;
	for( odc = 0 ; odc < NSNS ; odc++ ) {
	  // open dirac index for source and sink
	  const int OD1 = odc / NS ;
	  const int OD2 = odc % NS ;
	  // local accumulators
	  register double complex term1 = 0.0 ;
	  register double complex term2 = 0.0 ;
	  // Contract with the final propagator and trace out the source Dirac indices
	  // A polarization must still be picked for the two open Dirac indices offline
	  int dirac ;
	  for( dirac = 0 ; dirac < NS ; dirac++ ){
	    term1 += baryon_contract( DiQ , S2[ site ] , dirac , dirac , OD1 , OD2 ) ;
	    term2 += baryon_contract( DiQ , S2[ site ] , dirac , OD1 , dirac , OD2 ) ;
	  }
	  // Form the uds-, uud-type baryons and uuu-type distinguish from the Omega
	  Buds[ odc ] += term1 ;
	  Buud[ odc ] += term1 + term2 ;
	  Buuu[ odc ] += 2 * term1 + 4 * term2 ;
	}
      } // VOL3 loop

      // Fill baryon correlator array
      int i ;
      for( i = 0 ; i < NSNS ; i++ ) {
	Buds_corr[ GSRC ][ i ].C[ t ] = Buds[ i ] ;
	Buud_corr[ GSRC ][ i ].C[ t ] = Buud[ i ] ;
	Buuu_corr[ GSRC ][ i ].C[ t ] = Buuu[ i ] ;
      }
    }
    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_baryons( "Baryon: uds-type" , (const struct correlator**)Buds_corr ) ;
  debug_baryons( "Baryon: uud-type" , (const struct correlator**)Buud_corr ) ;
  debug_baryons( "Baryon: uuu-type" , (const struct correlator**)Buuu_corr ) ;
#endif

  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s.uds" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buds_corr , B_CHANNELS , NSNS ) ;
  sprintf( outstr , "%s.uud" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buud_corr , B_CHANNELS , NSNS ) ;
  sprintf( outstr , "%s.uuu" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buuu_corr , B_CHANNELS , NSNS ) ;

  free_corrs( Buds_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;

  // free the props
  free( S1 ) ;
  free( S2 ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}
