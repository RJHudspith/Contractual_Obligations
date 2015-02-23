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
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\rdone %.f %%",t/((L0-1)/100.));fflush(stdout);
  }
  printf("\n");	

  // & do something with the computed correlators
#ifdef DEBUG
  printf( "LL PION\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[5][5].C[t] ) , cimag( corr[5][5].C[t] ) ) ;
  }
  printf( "11\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[1][1].C[t] ) , cimag( corr[1][1].C[t] ) ) ;
  }
  printf( "1010\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[10][10].C[t] ) , cimag( corr[10][10].C[t] ) ) ;
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




// computes heavy-heavy meson correlators
// Note: Backwards propagators are automatically turned around due to the save-format
int
hheavy_mesons( FILE *prop1 , 
	       const int header )
{
  // data structure for holding the contractions
  struct correlator **corr = malloc( NS * NS * sizeof( struct correlator* ) ) ;
  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *P1 = malloc( VOL3 * sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_nrprop( prop1 , S1 , header , t ) ;

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
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\rdone %.f %%",t/((L0-1)/100.));fflush(stdout);
  }
  printf("\n");

  // & do something with the computed correlators
#ifdef DEBUG
  printf( "HH PION\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[5][5].C[t] ) , cimag( corr[5][5].C[t] ) ) ;
  }
  printf( "11\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[1][1].C[t] ) , cimag( corr[1][1].C[t] ) ) ;
  }
  printf( "1010\n" ) ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[10][10].C[t] ) , cimag( corr[10][10].C[t] ) ) ;
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
	       const int header,
	       const int header2 )
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
    read_nrprop( prop2 , S2 , header2 , t ) ;

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
          #ifdef DEBUG
	  sum += pion_correlator( S1[ site ] , S2[ site ] );
          #else 
	  sum += local_meson_correlator( S1[ site ] , S2[ site ] , 
					 GAMMAS[ GAMMA_FIVE ] , 
					 GAMMAS[ GSRC ] , 
					 GAMMAS[ GSNK ] ) ;
          #endif
	}
	//
	corr[ GAMMA_1 ][ GAMMA_2 ].C[ t ] = (double complex)sum ;
      }
    }
    printf("\rdone %.f %%",t/((L0-1)/100.));fflush(stdout);
  }
  printf("\n");

#ifdef DEBUG
  printf( "HL PION\n" ) ; 
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[5][5].C[t] ) , cimag( corr[5][5].C[t] ) ) ;
  } 
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
