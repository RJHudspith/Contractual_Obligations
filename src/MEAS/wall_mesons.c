/**
   @file wall_mesons.c
   @brief wall-wall meson coputation
 */

#include "common.h"

#include "correlators.h"  // for allocate_corrs and free_corrs
#include "gammas.h"       // gamma matrices
#include "io.h"           // read_prop
#include "spinor_ops.h"   // meson_contract

// sums a timeslice worth of S into SUM, is threaded but not how you would expect
static void
sumprop( struct spinor *__restrict SUM ,
	 struct spinor *__restrict S )
{
  int d1d2 ;
#pragma omp parallel for private(d1d2)
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {

    const int d1 = d1d2 / NS ;
    const int d2 = d1d2 % NS ;

    int c1 , c2 ;
    // initialise sum
    for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	SUM -> D[d1][d2].C[c1][c2] = 0.0 ;
      }
    }
    // memory access pattern is pretty strange
    int i ;
    for( i = 0 ; i < VOL3 ; i++ ) {
      // pull this into cache would do something like this ..
      // double complex *C = (double complex*)S[i].D[d1][d2].C ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  SUM -> D[d1][d2].C[c1][c2] += S[i].D[d1][d2].C[c1][c2] ;
	}
      }
    }
  }
  return ;
}

// computes meson correlators
int
wall_mesons( FILE *prop1 ,
	     const proptype proptype1 )
{
  printf( "Computing Wall-Source mesons \n" ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **wlcorr = malloc( NSNS * sizeof( struct correlator* ) ) ;
  struct correlator **wwcorr = malloc( NSNS * sizeof( struct correlator* ) ) ;

  allocate_corrs( wlcorr ) ;
  allocate_corrs( wwcorr ) ;

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 , proptype1 ) == FAILURE ) {
      free_corrs( wlcorr ) ;
      free_corrs( wwcorr ) ;
      free( GAMMAS ) ;
      free( S1 ) ;
      return FAILURE ;
    }
    
    // single sum
    struct spinor SUM1 ;
    sumprop( &SUM1 , S1 ) ;

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {
	
	register double complex sum = 0.0 ;

	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , 
				 S1[ site ] ,
				 GAMMAS[ GSRC ] ,
				 S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
				 
	}

	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;

	// correlator computed just out of the summed walls
	wwcorr[ GSRC ][ GSNK ].C[ t ] =		\
	  meson_contract( GAMMAS[ GSNK ] , 
			  SUM1 ,
			  GAMMAS[ GSRC ] ,
			  SUM1 ,
			  GAMMAS[ GAMMA_5 ] ) ;
      }
    }
  }

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
#endif

  // free our correlator measurement
  free_corrs( wlcorr ) ;

  // free our correlator measurement
  free_corrs( wwcorr ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;

  return SUCCESS ;
}

// computes meson correlators
int
wall_double_mesons( FILE *prop1 , 
		    const proptype proptype1 ,
		    FILE *prop2 ,
		    const proptype proptype2 )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **wlcorr = malloc( NSNS * sizeof( struct correlator* ) ) ;
  struct correlator **wwcorr = malloc( NSNS * sizeof( struct correlator* ) ) ;
  allocate_corrs( wlcorr ) ;
  allocate_corrs( wwcorr ) ;

  // and our spinors
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S2 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 , proptype1 ) == FAILURE ||
	read_prop( prop2 , S2 , proptype2 ) == FAILURE ) {
      free_corrs( wlcorr ) ;
      free_corrs( wwcorr ) ;
      free( GAMMAS ) ;
      free( S1 ) ;
      free( S2 ) ;
      return FAILURE ;
    }

    // prop sums
    struct spinor SUM1 , SUM2 ;
    sumprop( &SUM1 , S1 ) ;
    sumprop( &SUM2 , S2 ) ;

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {
	
	register double complex sum = 0.0 ;

	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , 
				 S2[ site ] ,
				 GAMMAS[ GSRC ] ,
				 S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}

	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;

	// correlator computed just out of the summed walls
	wwcorr[ GSRC ][ GSNK ].C[ t ] =	\
	  meson_contract( GAMMAS[ GSNK ] , 
			  SUM2 ,
			  GAMMAS[ GSRC ] ,
			  SUM1 ,
			  GAMMAS[ GAMMA_5 ] ) ;
      }
    }
  }

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
#endif

  // free our correlator measurement
  free_corrs( wlcorr ) ;

  // free our correlator measurement
  free_corrs( wwcorr ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  return SUCCESS ;
}
