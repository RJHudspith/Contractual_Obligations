/**
   @file mesonfile.c
   @brief little correlation file reader
 */
#include "common.h"

#include "correlators.h" // free_momcorr()
#include "reader.h"      // process_file()

// lattice geometry information
struct latt_info Latt ;

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 2 ) {
    return fprintf( stdout , "[MESONS] Usage ./MESONS {correlator file}"
		    " GSRC,GSNK,px,py,pz ... \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[MESONS] File %s does not exist\n" , argv[1] ) ;
    return FAILURE ;
  }

  // set geometry to 0
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 0 ;
  }

  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL ;
  struct veclist *momentum = NULL ;

  // number of correlators printed to stdout
  int corrs_written = 0 ;

  // read in the file and pack our structs
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;
  
  if( corr == NULL ) goto memfree ;

  // loop the ones we want
  size_t i ;
  for( i = 2 ; i < (size_t)( argc ) ; i++ ) {
    // tokenize argv into the correlators people want
    char *tok1 = strtok( (char*)argv[i] , "," ) ;
    if( tok1 == NULL ) break ;
    const int idx1 = atoi( tok1 ) ;
    if( idx1 >= (int)NGSRC[0] || idx1 < 0 ) { 
      fprintf( stderr , "[Momcorr] Non-sensical source index %d \n" , idx1 ) ;
      break ;
    } 
    char *tok2 = strtok( NULL , "," ) ;
    if( tok2 == NULL ) break ;
    const int idx2 = atoi( tok2 ) ;
    if( idx2 >= (int)NGSNK[0] || idx2 < 0 ) { 
      fprintf( stderr , "[Momcorr] Non-sensical sink index %d \n" , idx2 ) ;
      break ;
    } 

    // initialise to 0
    double moms[ ND - 1 ] ;
    int mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      moms[ mu ] = 0.0 ;
    }

    fprintf( stdout , "[Momcorr] searching for momentum (" ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      char *ptok = strtok( NULL , "," ) ;
      if( ptok == NULL ) break ;
      moms[ mu ] = (double)atof( ptok ) ;
      fprintf( stdout , " %g " , moms[ mu ] ) ;
    }
    fprintf( stdout , ") \n" ) ;

    // find the correlator in the list
    const size_t matchmom = find_desired_mom( momentum , moms , 
					      (int)NMOM[0] ) ;
    if( matchmom == 123456789 ) {
      fprintf( stderr , "[Momcorr] Unable to find desired momentum ... Leaving \n" ) ;
      break ;
    }

    fprintf( stdout , "[Momcorr] match ( %g %g %g ) \n" ,
	     momentum[ matchmom ].MOM[0] ,
	     momentum[ matchmom ].MOM[1] ,
	     momentum[ matchmom ].MOM[2] ) ;

    fprintf( stdout , "[Momcorr] Correlator [ Source :: %d | Sink :: %d ] \n\n" , 
	     idx1 , idx2 ) ;

    size_t t ;
    for( t = 0 ; t < LT ; t++ ) {
      fprintf( stdout , "CORR %zu %1.12e %1.12e\n" , t ,
	       creal( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ,
	       cimag( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ) ;
    }
    corrs_written++ ;
    //
    fprintf( stdout , "\n" ) ;
  }

  // if we don't have a match or didn't specify gammas give the momentum
  // list as an option
  if( corrs_written == 0 ) {
    write_momlist( momentum , NMOM[ 0 ] ) ;
  }

 memfree :

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  fclose( infile ) ;

  return SUCCESS ;
}
