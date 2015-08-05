/**
   @file baryonfile.c
   @brief little correlation file reader
   
   TODO :: momentum projections for spin 1/2, 3/2?
 */
#include "common.h"

#include "baryon_projections.h" // enums and projections
#include "correlators.h"        // free_momcorrs()
#include "gammas.h"             // make_gammas()
#include "reader.h"             // process file

// global lattice dimensions and stuff
struct latt_info Latt ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// little code for accessing elements of our baryon correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc != 6 ) {
    return printf( "usage ./BARYONS {correlator file} BASIS PROJECTION Lx,Ly,Lz {outfile} \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return FAILURE ;
  }
  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  struct gamma *GAMMAS = NULL ;     // Gamma matrices
  struct mcorr **corr = NULL ;      // read-in correlator
  struct veclist *momentum = NULL ; // momentum list
  struct mcorr **proj_corr = NULL ; // projected correlator

  // set the basis
  proptype basis = CHIRAL ;
  if( are_equal( argv[2] , "NREL" ) ) {
    basis = NREL ;
  } else if( are_equal( argv[2] , "STATIC" ) ) {
    basis = STATIC ;
  } else if( are_equal( argv[2] , "CHIRAL" ) ) {
    basis = CHIRAL ;
  } else {
    printf( "[INPUTS] I don't understand your basis %s \n" , argv[2] ) ;
    goto memfree ;
  }

  // set the projection
  bprojection projection = L0 ;
  if( are_equal( argv[3] , "L0" ) ) {
    printf( "[PARITY] Performing an L0 projection\n" ) ;
    projection = L0 ;
  } else if( are_equal( argv[3] , "L1" ) ) {
    printf( "[PARITY] Performing an L1 projection\n" ) ;
    projection = L1 ;
  } else if( are_equal( argv[3] , "L2" ) ) {
    printf( "[PARITY] Performing an L2 projection\n" ) ;
    projection = L2 ;
  } else if( are_equal( argv[3] , "L3" ) ) {
    printf( "[PARITY] Performing an L3 projection\n" ) ;
    projection = L3 ;
  } else if( are_equal( argv[3] , "L4" ) ) {
    printf( "[PARITY] Performing an L4 projection\n" ) ;
    projection = L4 ;
  } else if( are_equal( argv[3] , "L5" ) ) {
    printf( "[PARITY] Performing an L5 projection\n" ) ;
    projection = L5 ;
  } else {
    printf( "[INPUTS] I don't understand your projection %s \n" , argv[3] ) ;
    goto memfree ;
  }

  // get the dimensions
  int mu ;
  char *tok = strtok( (char*)argv[4] , "," ) ;
  Latt.dims[ 0 ] = (int)atoi( tok ) ;
  if( Latt.dims[ 0 ] < 1 ) {
    printf( "[INPUTS] non-sensical lattice dimension %d %d \n" , 
	    0 , Latt.dims[0] ) ;
    goto memfree ;
  }
  for( mu = 1 ; mu < ND-1 ; mu++ ) {
    char *ptok = strtok( NULL , "," ) ;
    if( ptok == NULL ) break ;
    Latt.dims[ mu ] = (int)atoi( ptok ) ;
    if( Latt.dims[mu] < 1 ) {
      printf( "[INPUTS] non-sensical lattice dimension %d %d \n" , 
	      mu , Latt.dims[mu] ) ;
      goto memfree ;
    }
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , basis ) == FAILURE ) {
    goto memfree ;
  }

  // is defined in reader.c
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;

  if( corr == NULL ) goto memfree ;

  // allocate corr
  
  proj_corr = allocate_momcorrs( B_CHANNELS , B_CHANNELS , NMOM[0] ) ;

  // do the projection
  size_t GSGK ;
#pragma omp parallel for private( GSGK )
  for( GSGK = 0 ; GSGK < ( B_CHANNELS * B_CHANNELS ) ; GSGK++ ) {

    const size_t GSRC = GSGK / B_CHANNELS ;
    const size_t GSNK = GSGK % B_CHANNELS ;

    // loop momenta
    size_t p ;
    for( p = 0 ; p < NMOM[0] ; p++ ) {

      // projections happen here
      const double complex *C = baryon_project( (const struct mcorr**)corr , 
						GAMMAS , momentum ,
						GSRC , GSNK , p ,
						projection ) ;
      
      // poke into proj_corr
      size_t t ;
      for( t = 0 ; t < LT ; t++ ) {
	proj_corr[ GSRC ][ GSNK ].mom[ p ].C[ t ] = C[ t ] ;
      }
    
      free( (void*)C ) ;
    }
  }

  // write out the correlator
  write_momcorr( argv[5] , (const struct mcorr **)proj_corr , 
		 momentum , B_CHANNELS , B_CHANNELS , (const int*)NMOM ) ;

 memfree :

  // free the gamma matrices
  free( GAMMAS ) ;

  // free the memory of the read-in correlator
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the memory of the projected correlator
  free_momcorrs( proj_corr , B_CHANNELS , B_CHANNELS , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  // close the infile
  fclose( infile ) ;

  return 0 ;
}
