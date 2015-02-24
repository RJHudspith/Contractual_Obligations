/**
   @file read_propheader.c
   @brief propagator header reader
 */
#include "common.h"

static int
are_equal( const char *pch , const char *tag )
{
  return !strcmp( pch , tag ) ;
}

static char*
get_tag( char line[ MAX_LINE_LENGTH ] )
{
  return strtok( line , " " ) ;
}

//
static int
read_header( FILE *prop ,
	     int dims[ ND ] ,
	     int src[ ND ] )
{
  const int MAX_HEADER_LINES = 32 ;
  char line[ MAX_LINE_LENGTH ] ;
  int n = 0 ;
  while( n < MAX_HEADER_LINES ) {
    if( fgets( line , MAX_LINE_LENGTH , prop ) == NULL ) {
      printf( "[IO] prop Header reading failed ... Leaving \n" ) ;
      return FAILURE ;
    }
    char *tag = get_tag( line ) ;

    // tag search for Lattice
    if( are_equal( tag , "Lattice:" ) ) {
      int N = 0 ;
      char *token ;
      while( ( token = strtok( NULL , " " ) ) != NULL ) {
	dims[ N ] = (int)atoi( token ) ;
	if( N == ND ) break ;
	N++ ;
      }
    }

    // tag search for source position
    if( are_equal( tag , "SrcPos:" ) ) {
      int N = 0 ;
      char *token ;
      while( ( token = strtok( NULL , " " ) ) != NULL ) {
	src[ N ] = (int)atoi( token ) ;
	if( N == ND ) break ;
	N++ ;
      }
    }

    // could put the other stuff here in good time

    if( strcmp ( line , "<end_header>\n" ) == 0 ) {
      break ;
    }
    n++ ;
  }
  return SUCCESS ;
}

int
read_check_header( FILE *propfile )
{
  int dims[ ND ] = {} ;
  int src[ ND ] ;

  if( read_header( propfile , dims , src ) == FAILURE ) {
    return FAILURE ;
  }

  // check lattice dimensions
  printf( "PROP DIMS :: ( " ) ;
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( dims[ mu ] != Latt.dims[ mu ] ) {
      printf( "prop header dimensions do not match lattice dims .. Leaving\n" ) ;
      return FAILURE ;
    }
    printf( "%d " , dims[ mu ] ) ;
  }
  printf( ")\n" ) ;

  // maybe we should check the source position as well?
  
  return SUCCESS ;
}
