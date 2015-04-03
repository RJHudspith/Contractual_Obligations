/**
   @file baryons.c
   @brief baryon contraction code


   This needs some splainin
*/

#include "common.h"

#include "bar_contractions.h"  // cross_color trace and baryon contract
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // allocate_corrs() && free_corrs()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas() && gamma_mmul
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()

// flavour degenerate baryon contraction
int
baryons_diagonal( struct propagator prop ,
		  const char *outfile )
{
  // Define our output correlators, with B_CHANNELS channels and NSNS components
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

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free_corrs( Buds_corr , 6 , NSNS ) ; 
    free_corrs( Buud_corr , 6 , NSNS ) ; 
    free_corrs( Buuu_corr , 6 , NSNS ) ; 
    free( S1 ) ; free( GAMMAS ) ;
    return FAILURE ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop , S1 ) == FAILURE ) {
      free_corrs( Buds_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buud_corr , B_CHANNELS , NSNS ) ; 
      free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ; 
      free( S1 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop.basis == CHIRAL && prop.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } 

    int GSRC = 0 ;
    // parallelise the furthest out loop
#pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {

      // Define some intermediate spinors
      struct spinor DiQ ;

      // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
      const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
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

	// Cross color product and sink Dirac trace
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
	    term1 += baryon_contract( DiQ, S1[ site ] , dirac , dirac , OD1 , OD2 ) ;
	    term2 += baryon_contract( DiQ, S1[ site ] , dirac , OD1 , dirac , OD2 ) ;
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

  free( S1 ) ;
  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}

#if 0
/*
  struct spinor CgS, CgS51, CgS52 ;	
  struct spinor DiQ, DiQ51, DiQ52 ;	

  term3 = 0.0 ;
  term4 = 0.0 ;
  term5 = 0.0 ;
  term6 = 0.0 ;

  Omega stuff
  // In case of the Omega ( s ( s Cg5 s )) we need the source and sink separately
  CgS51 = Cgamma_snk( S1[ site ] , GSRC ); -> gamma_mul_l( S1[ site ] , CgmuT ) 
  CgS52 = Cgamma_src( S1[ site ] , GSRC ); -> gamma_mul_r( S1[ site ] , Cgmu ) 

  // There are four more terms for the Omega					
  term3 += baryon_contract( DiQ51, S1[ site ], dirac , dirac , OD1 , OD2 );
  term4 += baryon_contract( DiQ52, CgS52, dirac , OD1 , dirac , OD2 );
  term5 += baryon_contract( DiQ52, CgS52, OD1 , dirac , dirac , OD2 );
  term6 += baryon_contract( DiQ51, S1[ site ], OD1 , dirac , dirac , OD2 );

  Buuu += term1 + term2 + term3 + term4 + term5 + term6;

  // In case of the Omega we need two additional diquarks
  DiQ51 = cross_color_trace( CgS52, CgS51 );
  DiQ52 = cross_color_trace( S1[ site ], CgS51 );
 */
#endif
