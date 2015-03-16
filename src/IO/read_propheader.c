/**
   @file read_propheader.c
   @brief propagator header reader
 */
#include "common.h"

#include <errno.h>    // for the error codes to strto*

//#define NEWHEADER

// dirty string equivalence
static int
are_equal( const char *pch , const char *tag )
{
  return !strcmp( pch , tag ) ;
}

// get a particular tag
static char*
get_tag( char line[ MAX_LINE_LENGTH ] )
{
  return strtok( line , " " ) ;
}

// get the propagator dims
static int
get_propdims( int *dims ) 
{
  char *endptr , *token ;
  int N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    dims[ N ] = (int)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE || dims[ N ] < 0 ) {
      printf( "[IO] propheader dims[%d] misread %d \n" , N , dims[ N ] ) ;
      return FAILURE ;
    }
    // check against global lattice data struct
    if( dims[ N ] != Latt.dims[ N ] ) {
      printf( "[IO] propheader and global lattice dims mismatch!\n"
	      "[IO] %d vs %d ( index %d ) \n" , dims[ N ] ,
	      Latt.dims[ N ] , N ) ;
    }
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
} 

// get the initial source position
static int
get_propsrc( int *origin ) 
{
  char *endptr , *token ;
  int N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    origin[ N ] = (int)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE || 
	origin[ N ] < 0 || origin[ N ] > Latt.dims[ N ] ) {
      printf( "[IO] propheader SrcPos:[%d] misread %d \n" , N , origin[ N ] ) ;
      return FAILURE ;
    }
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

// get the propagator endianness
static int
get_propendian( endianness *endian )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Big" ) ) {
      *endian = BIGENDIAN ;
      return SUCCESS ;
    }
    if( are_equal( token , "Little" ) ) {
      *endian = LILENDIAN ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// get the propagator precision
static int
get_propprec( fp_precision *precision )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Single" ) ) {
      *precision = SINGLE ;
      return SUCCESS ;
    }
    if( are_equal( token , "Double" ) ) {
      *precision = DOUBLE ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// get the propagator source type
static int
get_propsource( sourcetype *source )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Point" ) ) {
      *source = POINT ;
      return SUCCESS ;
    }
    // wall can be a bunch of things I think
    if( are_equal( token , "GFWall" ) ||
	are_equal( token , "Z2Wall" ) ||
	are_equal( token , "Wall" ) ) {
      *source = WALL ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// get the propagator source type
static int
get_proptype( proptype *basis )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // nonrelativistic basis, this may need to become
    // Nrel_fwd and Nrel_bwd 
    if( are_equal( token , "Nrel" ) ) {
      *basis = NREL ;
      return SUCCESS ;
    }
    // relativistic chiral basis
    if( are_equal( token , "Chiral" ) ) {
      *basis = CHIRAL ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// little sugar
static int
tagfailure( const char *message , 
	    const char *line )
{
  printf( "[IO] propheader failure while searching for tag %s \n"
	  "[IO] line in file is %s \n" , message , line ) ;
  return FAILURE ;
}

// tells us if we missed a record we were looking for
static int
nonexistent_record( const char *message )
{
  printf( "[IO] I cannot find record for tag %s \n" , message ) ;
  return FAILURE ;
}

//
int
read_propheader( struct propagator *prop )
{
  // dimensions in the file
  int dims[ ND ] = { } ;
  const int MAX_HEADER_LINES = 32 ;
  char line[ MAX_LINE_LENGTH ] , *endptr ;
  int n = 0 , dimsflag = 0 , originflag = 0 ;
#ifdef NEWHEADER
  int endflag = 0 , srcflag = 0 , precflag = 0 ;
  int basisflag = 0 ;
#endif

  // loop up to some large number of possible tags
  while( n < MAX_HEADER_LINES ) {

    // read a line
    if( fgets( line , MAX_LINE_LENGTH , prop -> file ) == NULL ) {
      printf( "[IO] prop Header reading failed ... Leaving \n" ) ;
      return FAILURE ;
    }
    char *tag = get_tag( line ) ;

    // tag search for Lattice
    if( are_equal( tag , "Lattice:" ) ) {
      if( get_propdims( dims ) == FAILURE ) {
	return tagfailure( "Lattice:" , line ) ;
      }
      dimsflag ++ ; 
    }
    // tag search for source origin
    if( are_equal( tag , "SrcPos:" ) ) {
      if( get_propsrc( prop -> origin ) == FAILURE ) {
	return tagfailure( "SrcPos:" , line ) ;
      }
      originflag ++ ; 
    }
#ifdef NEWHEADER
    // incorporate these when the header catches up
    if( are_equal( tag , "Endian:" ) ) {
      if( get_propendian( &( prop -> endian ) ) == FAILURE ) {
	return tagfailure( "Endian:" , line ) ;
      }
      endflag ++ ;
    }
    // precision
    if( are_equal( tag , "Precision:" ) ) {
      if( get_propprec( &( prop -> precision ) ) == FAILURE ) {
	return tagfailure( "Precision:" , line ) ;
      }
      precflag ++ ;
    }
    // source
    if( are_equal( tag , "Source:" ) ) {
      if( get_propsource( &( prop -> source ) ) == FAILURE ) {
	return tagfailure( "Source:" , line ) ;
      }
      srcflag ++ ;
    }
    // basis
    if( are_equal( tag , "Basis:" ) ) {
      if( get_proptype( &( prop -> basis ) ) == FAILURE ) {
	return tagfailure( "Basis:" , line ) ;
      }
      basisflag ++ ;
    }
#endif

    if( are_equal( line , "<end_header>\n" ) ) {
      break ;
    }
    n++ ;
  }

  // check the flags have been hit
  if( dimsflag == 0 ) return nonexistent_record( "Lattice:" ) ; 
  if( originflag == 0 ) return nonexistent_record( "SrcPos:" ) ; 
#ifdef NEWHEADER
  if( endflag == 0 ) return nonexistent_record( "Endian:" ) ; 
  if( precflag == 0 ) return nonexistent_record( "Precision:" ) ;
  if( srcflag == 0 ) return nonexistent_record( "Source:" ) ;
  if( basisflag == 0 ) return nonexistent_record( "Basis:" ) ;
#endif
  if( n == MAX_HEADER_LINES ) return nonexistent_record( "<end header>" ) ;

  // if everything is kosher we leave successfully
  return SUCCESS ;
}

// open the files and parse the header
int
read_propheaders( struct propagator *prop ,
		  const struct input_info inputs )
{
  int i ;
  for( i = 0 ; i < inputs.nprops ; i++ ) {
    // read and check 'em
    if( read_propheader( &prop[ i ] ) == FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
