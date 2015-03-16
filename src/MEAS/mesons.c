/**
   @file wall_mesons.c
   @brief wall-wall meson coputation
 */

#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time() 
#include "io.h"                // read_prop
#include "spinor_ops.h"        // spinor operations
#include "read_propheader.h"   // (re)read the propagator header

// computes meson correlators
int
mesons_diagonal( struct propagator prop ,
		 const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **wlcorr = allocate_corrs( NSNS , NSNS ) ;
  struct correlator **wwcorr = NULL ;

  if( prop.source == WALL ) {
    wwcorr = allocate_corrs( NSNS , NSNS ) ;
  }

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop , S1 ) == FAILURE ) {
      free_corrs( wlcorr , NSNS , NSNS ) ; 
      if( prop.source == WALL ) {
	free_corrs( wwcorr , NSNS , NSNS ) ;
      }
      free( GAMMAS ) ; free( S1 ) ;
      return FAILURE ;
    }
    
    // single sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {

	register double complex sum = 0.0 ;

	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S1[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}

	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;

	if( prop.source == WALL ) {
	  // correlator computed just out of the summed walls
	  wwcorr[ GSRC ][ GSNK ].C[ t ] =		\
	    meson_contract( GAMMAS[ GSNK ] , SUM1 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
	//
      }
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  if( prop.source == WALL ) {
    debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
  }
#endif

  // and write the output files
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)wlcorr ,
		     NSNS , NSNS ) ;

  // free our correlator measurement
  free_corrs( wlcorr , NSNS , NSNS ) ;

  // if we have computed the wall-wall we output it
  if( prop.source == WALL ) {
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)wwcorr ,
		       NSNS , NSNS ) ;
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}

// computes meson correlators
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **wlcorr = allocate_corrs( NSNS , NSNS ) ;
  struct correlator **wwcorr = NULL ;
  
  if( prop1.source == WALL || prop2.source == WALL ) {
    allocate_corrs( NSNS , NSNS ) ;
  }

  // and our spinors
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S2 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 ) == FAILURE ||
	read_prop( prop2 , S2 ) == FAILURE ) {
      free_corrs( wlcorr , NSNS , NSNS ) ; 
      if( prop1.source == WALL || prop2.source == WALL ) {
	free_corrs( wwcorr , NSNS , NSNS ) ;
      }
      free( GAMMAS ) ; free( S1 ) ; free( S2 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
      nrel_rotate_slice( S2 ) ;
    }

    // prop sums
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {

	register double complex sum = 0.0 ;

	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S2[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}

	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;

	if( prop1.source == WALL || prop2.source == WALL ) {
	  // correlator computed just out of the summed walls
	  wwcorr[ GSRC ][ GSNK ].C[ t ] =		\
	    meson_contract( GAMMAS[ GSNK ] , SUM2 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
	//
      }
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  if( prop1.source == WALL || prop2.source == WALL ) {
    debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
  }
#endif

  // outputs
  char outstr[ 256 ] ;
  sprintf( outstr , "%s.wl" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)wlcorr ,
		     NSNS , NSNS ) ;

  free_corrs( wlcorr , NSNS , NSNS ) ;

  if( prop1.source == WALL || prop2.source == WALL ) {
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)wwcorr ,
		       NSNS , NSNS ) ;
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}
