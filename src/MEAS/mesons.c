/**
   @file mesons.c
   @brief meson computation codes

   TODO :: write out computed correlated matrix, format? crc?
 */

#include "common.h"

#include "correlators.h"  // correlators
#include "gammas.h"       // gamma matrices
#include "io.h"           // read prop

#define DEBUG

// allocate the correlator matrix
static void
allocate_corrs( struct correlator **corr )
{
  int s ;
  for( s = 0 ; s < NS * NS ; s++ ) {
    corr[ s ] = ( struct correlator* )malloc( NS * NS * sizeof( struct correlator ) ) ;
    int t ;
    for( t = 0 ; t < NS * NS ; t++ ) {
      corr[s][t].C = malloc( L0 * sizeof( double complex ) ) ;
    }
  }
  return ;
}

// free the correlator matrix
static void
free_corrs( struct correlator **corr )
{
  int s ;
  for( s = 0 ; s < NS * NS ; s++ ) {
    int t ;
    for( t = 0 ; t < NS * NS ; t++ ) {
      free( corr[ s ][ t ].C ) ;
    }
    free( corr[ s ] ) ;
  }
  free( corr ) ;
  return ;
}

// computes meson correlators
int
single_mesons( FILE *prop1 , 
	       const int header )
{
  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;

  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , header , t ) ;

    int G1 = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(G1)
    for( G1 = 0 ; G1 < NS*NS ; G1++ ) {

      int G2 ;
      for( G2 = 0 ; G2 < NS*NS ; G2++ ) {
	
	register double complex sum = 0.0 ;

	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += local_meson_correlator_singlet( S1[ site ] , 
						 GAMMAS[ GAMMA_5 ] , 
						 GAMMAS[ G1 ] , 
						 GAMMAS[ G2 ] ) ;
	}
	//
	corr[ G1 ][ G2 ].C[ t ] = (double complex)sum ;
      }
    }
  }

  // & do something with the computed correlators
#ifdef DEBUG
  printf( "PION\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_5][GAMMA_5].C[t] ) , 
	    cimag( corr[GAMMA_5][GAMMA_5].C[t] ) ) ;
  }
  printf( "00\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_0][GAMMA_0].C[t] ) , 
	    cimag( corr[GAMMA_0][GAMMA_0].C[t] ) ) ;
  }
  printf( "11\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_1][GAMMA_1].C[t] ) , 
	    cimag( corr[GAMMA_1][GAMMA_1].C[t] ) ) ;
  }
  printf( "22\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_2][GAMMA_2].C[t] ) , 
	    cimag( corr[GAMMA_2][GAMMA_2].C[t] ) ) ;
  }
  printf( "33\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_3][GAMMA_3].C[t] ) , 
	    cimag( corr[GAMMA_3][GAMMA_3].C[t] ) ) ;
  }
#endif

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;

  return SUCCESS ;
}

// computes meson correlators from two propagators
int
double_mesons( FILE *prop1 , 
	       FILE *prop2 ,
	       const int header )
{
  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;
 
  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S2 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , header , t ) ;
    read_prop( prop2 , S2 , header , t ) ;

    int G1 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(G1)
    for( G1 = 0 ; G1 < NS*NS ; G1++ ) {

      int G2 ;
      for( G2 = 0 ; G2 < NS*NS ; G2++ ) {

	register complex sum = 0.0 ;
	//
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += local_meson_correlator( S1[ site ] , S2[ site ] , 
					 GAMMAS[ GAMMA_5 ] , 
					 GAMMAS[ G1 ] , 
					 GAMMAS[ G2 ] ) ;
	}
	//
	corr[ G1 ][ G2 ].C[ t ] = (double complex)sum ;
      }
    }
  }

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  return SUCCESS ;
}
