/**
   @file mesons.c
   @brief meson computation codes

   TODO :: write out computed correlated matrix, format? crc?
 */

#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print time
#include "io.h"                // read prop
#include "read_propheader.h"   // (re)read the header

// computes meson correlators
int
single_mesons( FILE *prop1 , 
	       const proptype proptype1 ,
	       const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;

  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 , proptype1 ) == FAILURE ) {
      free_corrs( corr ) ;
      free( GAMMAS ) ;
      free( S1 ) ;
      return FAILURE ;
    }

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NS*NS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NS*NS ; GSNK++ ) {

	register double complex sum = 0.0 ;

	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S1[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	//
	corr[ GSRC ][ GSNK ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\r[MESONS] done %.f %%",(t+1)/((L0)/100.));fflush(stdout);
  }
  printf("\n");	

  // & do something with the computed correlators
#ifdef DEBUG
  debug_mesons( "LL-mesons" , (const struct correlator**)corr ) ;
#endif

  // and write out a file
  write_correlators( outfile , (const struct correlator**)corr ) ;

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;

  // rewind file and read header again
  rewind( prop1 ) ;
  read_check_header( prop1 , GLU_FALSE ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}

// computes meson correlators from two propagators
int
double_mesons( FILE *prop1 , 
	       const proptype proptype1 ,
	       FILE *prop2 ,
	       const proptype proptype2 ,
	       const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the non relativistic gamma basis
  if( make_gammas( GAMMAS , proptype1 ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // data structure for holding the contractions
  struct correlator **corr = malloc( NSNS * sizeof( struct correlator* ) ) ;
 
  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S2 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop1 , S1 , proptype1 ) == FAILURE ||
	read_prop( prop2 , S2 , proptype2 ) == FAILURE ) {
      free_corrs( corr ) ;
      free( GAMMAS ) ;
      free( S1 ) ;
      free( S2 ) ;
      return FAILURE ;
    }
    
    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( proptype1 == CHIRAL && proptype2 == NREL ) {
      int site ;
      #pragma omp parallel for private(site) 
      for( site = 0 ; site < LCU ; site++ ) {
	chiral_to_nrel( &S1[ site ] ) ;
      }
    } else if( proptype1 == NREL && proptype2 == CHIRAL ) {
      int site ;
      #pragma omp parallel for private(site) 
      for( site = 0 ; site < LCU ; site++ ) {
	chiral_to_nrel( &S2[ site ] ) ;
      }
    }

    int GSRC ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NS*NS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NS*NS ; GSNK++ ) {

	register double complex sum = 0.0 ;

	//
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S2[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	//
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; fflush( stdout ) ;
  }
  printf("\n");

#ifdef DEBUG
  debug_mesons( "HL-mesons" , (const struct correlator**)corr ) ;
#endif

  // and write out a file
  write_correlators( outfile , (const struct correlator**)corr ) ;

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  // rewind file and read header again
  rewind( prop1 ) ;
  read_check_header( prop1 , GLU_FALSE ) ;
  rewind( prop2 ) ;
  read_check_header( prop2 , GLU_FALSE ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}
