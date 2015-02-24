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
void
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
void
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

#ifdef DEBUG
void
debug_mesons( const char *message , 
	      const struct correlator **corr )
{
  int t ;
  printf( "%s PION\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_5][GAMMA_5].C[t] ) , cimag( corr[GAMMA_5][GAMMA_5].C[t] ) ) ;
  }
  printf( "%s 00\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_0][GAMMA_0].C[t] ) , cimag( corr[GAMMA_0][GAMMA_0].C[t] ) ) ;
  }
  printf( "%s 11\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_1][GAMMA_1].C[t] ) , cimag( corr[GAMMA_1][GAMMA_1].C[t] ) ) ;
  }
  printf( "%s 22\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_2][GAMMA_2].C[t] ) , cimag( corr[GAMMA_2][GAMMA_2].C[t] ) ) ;
  }
  printf( "%s 33\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GAMMA_3][GAMMA_3].C[t] ) , cimag( corr[GAMMA_3][GAMMA_3].C[t] ) ) ;
  }
  printf( "%s 1010\n" , message ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[10][10].C[t] ) , cimag( corr[10][10].C[t] ) ) ;
  }
}
#endif

// computes meson correlators
int
single_mesons( FILE *prop1 , 
	       const proptype proptype1 )
{
  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;

  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS , proptype1 ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , proptype1 ) ;

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
	  sum += local_meson_correlator( S1[ site ] , S1[ site ] , 
					 GAMMAS[ GAMMA_5 ] , 
					 GAMMAS[ GSRC ] , 
					 GAMMAS[ GSNK ] ) ;
	}
	//
	corr[ GSRC ][ GSNK ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\rdone %.f %%",t/((L0-1)/100.));fflush(stdout);
  }
  printf("\n");	

  // & do something with the computed correlators
#ifdef DEBUG
  debug_mesons( "LL-mesons" , (const struct correlator**)corr ) ;
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
	       const proptype proptype1 ,
	       FILE *prop2 ,
	       const proptype proptype2 )
{
  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;
 
  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S2 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the non relativistic gamma basis
  make_gammas( GAMMAS , proptype1 ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , proptype1 ) ;
    read_prop( prop2 , S2 , proptype2 ) ;

    int GSRC ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NS*NS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NS*NS ; GSNK++ ) {

	register complex sum = 0.0 ;
	//
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += local_meson_correlator( S1[ site ] , S2[ site ] , 
					 GAMMAS[ GAMMA_5 ] , 
					 GAMMAS[ GSRC ] , 
					 GAMMAS[ GSNK ] ) ;
	}
	//
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\rdone %.f %%",t/((L0-1)/100.));fflush(stdout);
  }
  printf("\n");

#ifdef DEBUG
  debug_mesons( "HL-mesons" , (const struct correlator**)corr ) ;
#endif

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  return SUCCESS ;
}
