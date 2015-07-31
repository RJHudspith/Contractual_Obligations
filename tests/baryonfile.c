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
  if( argc < 4 ) {
    return printf( "usage ./BARYONS {correlator file} BASIS PROJECTION Lx,Ly,Lz GSRC,GSNK,px,py,pz ... \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return FAILURE ;
  }
  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // structs
  struct gamma *GAMMAS = NULL ;
  struct mcorr **corr = NULL ;
  struct veclist *momentum = NULL ;

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
    projection = L0 ;
  } else if( are_equal( argv[3] , "L1" ) ) {
    projection = L1 ;
  } else if( are_equal( argv[3] , "L2" ) ) {
    projection = L2 ;
  } else if( are_equal( argv[3] , "L3" ) ) {
    projection = L3 ;
  } else if( are_equal( argv[3] , "L4" ) ) {
    projection = L4 ;
  } else if( are_equal( argv[3] , "L5" ) ) {
    projection = L5 ;
  } else {
    printf( "[INPUTS] I don't understand your projection %s \n" , argv[3] ) ;
    goto memfree ;
  }

  // get the dimensions
  int mu ;
  char *tok = strtok( (char*)argv[4] , "," ) ;
  Latt.dims[ 0 ] = (int)atoi( tok ) ;
  for( mu = 1 ; mu < ND-1 ; mu++ ) {
    char *ptok = strtok( NULL , "," ) ;
    if( ptok == NULL ) break ;
    Latt.dims[ mu ] = (int)atoi( ptok ) ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , basis ) == FAILURE ) {
    goto memfree ;
  }

  // number of correlators printed to stdout
  int corrs_written = 0 ;

  // is defined in reader.c
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;

  if( corr == NULL ) goto memfree ;

  // loop the ones we want
  int i ;
  for( i = 5 ; i < ( argc ) ; i++ ) {
    // tokenize argv into the correlators people want
    char *tok1 = strtok( (char*)argv[i] , "," ) ;
    if( tok1 == NULL ) break ;
    const int idx1 = atoi( tok1 ) ;
    if( idx1 >= B_CHANNELS || idx1 < 0 ) { 
      printf( "[Momcorr] Non-sensical source index %d \n" , idx1 ) ;
      break ;
    } 
    char *tok2 = strtok( NULL , "," ) ;
    if( tok2 == NULL ) break ;
    const int idx2 = atoi( tok2 ) ;
    if( idx2 >= B_CHANNELS || idx2 < 0 ) { 
      printf( "[Momcorr] Non-sensical sink index %d \n" , idx2 ) ;
      break ;
    } 

    // initialise to 0
    int moms[ ND - 1 ] , mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      moms[ mu ] = 0 ;
    }

    printf( "[Momcorr] searching for momenta (" ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      char *ptok = strtok( NULL , "," ) ;
      if( ptok == NULL ) break ;
      moms[ mu ] = (int)atoi( ptok ) ;
      printf( " %d " , moms[ mu ] ) ;
    }
    printf( ") \n" ) ;

    // find the correlator in the list
    const int matchmom = find_desired_mom( momentum , moms , 
					   (int)NMOM[0] ) ;
    if( matchmom == FAILURE ) {
      printf( "[Momcorr] Unable to find desired momentum ... Leaving \n" ) ;
      break ;
    }

    printf( "[Momcorr] match ( %d %d %d ) \n" , momentum[ matchmom ].MOM[ 0 ] ,
	    momentum[ matchmom ].MOM[ 1 ] ,  momentum[ matchmom ].MOM[ 2 ] ) ;

    printf( "[Momcorr] Correlator [ Source :: %d | Sink :: %d ] \n\n" , 
	    idx1 , idx2 ) ;

    // do the projection
    const double complex *C = baryon_project( (const struct mcorr**)corr , 
					      GAMMAS , 
					      momentum ,
					      idx1 , idx2 , matchmom ,
					      projection ) ;
    
    size_t t ;
    for( t = 0 ; t < LT ; t++ ) {
      printf( "CORR %zu %1.12e %1.12e\n" , t ,
	      creal( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ,
	      cimag( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ) ;
    }

    free( (void*)C ) ;
    corrs_written++ ;
    //
    printf( "\n" ) ;
  }

  // if we don't have a match or didn't specify gammas give the momentum
  // list as an option
  if( corrs_written == 0 ) {
    write_momlist( momentum , NMOM[ 0 ] ) ;
  }

 memfree :

  // free the gamma matrices
  free( GAMMAS ) ;

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  fclose( infile ) ;

  return 0 ;
}
