/**
   @file mesons.c
   @brief meson computation codes
 */

#include "common.h"

#include "correlators.h"  // correlators
#include "gammas.h"       // gamma matrices
#include "io.h"           // read prop

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

  // Compute lookup index for adjoint props 
  int iadj[ NS ] , sadj[ NS ] ;
  gamma_matrix( sadj , iadj , 5 ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , header , t ) ;

    int GAMMA_1 = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GAMMA_1)
    for( GAMMA_1 = 0 ; GAMMA_1 < NS*NS ; GAMMA_1++ ) {
      int GAMMA_2 ;
      for( GAMMA_2 = 0 ; GAMMA_2 < NS*NS ; GAMMA_2++ ) {
	register double complex sum = 0.0 ;
	//
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += local_meson_correlator( S1[ site ] , S1[ site ] , 
					 sadj , iadj , 
					 GAMMA_1 , GAMMA_2 ) ;
	}
	//
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
  }

  // free our correlator measurement
  free_corrs( corr ) ;

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

  // Compute lookup index for adjoint props 
  int iadj[ NS ] , sadj[ NS ] ;
  gamma_matrix( sadj , iadj , 5 ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , header , t ) ;
    read_prop( prop2 , S2 , header , t ) ;

    int GAMMA_1 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GAMMA_1)
    for( GAMMA_1 = 0 ; GAMMA_1 < NS*NS ; GAMMA_1++ ) {
      int GAMMA_2 ;
      //for( GAMMA_2 = 0 ; GAMMA_2 < NS*NS ; GAMMA_2++ ) {
      for( GAMMA_2 = 0 ; GAMMA_2 < 1 ; GAMMA_2++ ) {

	register complex sum = 0.0 ;
	//
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum = sum + (double complex)local_meson_correlator( S1[ site ] , S2[ site ] , 
							      sadj , iadj , 
							      GAMMA_1 , GAMMA_2 ) ;
	}
	//
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
  }

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  return SUCCESS ;
}
