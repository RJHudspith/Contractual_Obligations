/**
   @file read_propheader.c
   @brief propagator header reader

   The way this works is the code reads a line associated with a "tag"

   e.g.
   Lattice: {stuff}

   If the tag is something we want it calls a get_* function
   these then increment whitespace until the end of the line looking for
   what we wish in the case of the lattice it will look for data of the form

   DIM_X DIM_Y DIM_Z DIM_T and compare to the global geometry

   In the case of places where a single value is expected it returns
   as soon as it has found the value

   e.g.
   Endian: Big Little Medium
   Would return Big
 */
#include "common.h"

#include <errno.h>    // for the error codes to strto*

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
get_propdims( size_t *dims ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    dims[ N ] = (size_t)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] propheader dims[%zu] misread %zu \n" , 
	       N , dims[ N ] ) ;
      return FAILURE ;
    }
    // check against global lattice data struct
    if( dims[ N ] != Latt.dims[ N ] ) {
      fprintf( stderr , "[IO] propheader and global lattice dims mismatch!\n"
	       "[IO] %zu vs %zu ( index %zu ) \n" , dims[ N ] ,
	       Latt.dims[ N ] , N ) ;
      return FAILURE ;
    }
    // we have successfully read in all the dimension info we need
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
} 

// get the initial source position
static int
get_propsrc( size_t *origin ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    origin[ N ] = (size_t)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE || 
	origin[ N ] > Latt.dims[ N ] ) {
      fprintf( stderr , "[IO] propheader SrcPos:[%zu] misread %zu \n" , 
	       N , origin[ N ] ) ;
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
    if( are_equal( token , "Big\n" ) ) {
      *endian = BIGENDIAN ;
      return SUCCESS ;
    }
    if( are_equal( token , "Little\n" ) ) {
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
    if( are_equal( token , "Single\n" ) ) {
      *precision = SINGLE ;
      return SUCCESS ;
    }
    if( are_equal( token , "Double\n" ) ) {
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
    if( are_equal( token , "Point\n" ) || 
	are_equal( token , "Pt\n" ) || 
	are_equal( token , "pt\n" ) ) {
      *source = POINT ;
      return SUCCESS ;
    }
    // wall can be a bunch of things I think
    if( are_equal( token , "GFWall\n" ) ||
	are_equal( token , "Z2Wall\n" ) ||
	are_equal( token , "Wall\n" ) ||
	are_equal( token , "wall\n" ) || 
	are_equal( token , "CGwall\n" ) ) {
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
    if( are_equal( token , "Nrel_fwd\n" ) ) {
      *basis = NREL_FWD ;
      return SUCCESS ;
    } else if( are_equal( token , "Nrel_bwd\n" ) ) {
      *basis = NREL_BWD ;
      return SUCCESS ;
    } else if( are_equal( token , "Static\n" ) ) {
      *basis = STATIC ;
      return SUCCESS ;
    } else if( are_equal( token , "Chiral\n" ) ) {
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
  fprintf( stderr , "[IO] propheader failure while searching for tag %s \n"
	   "[IO] line in file is %s \n" , message , line ) ;
  return FAILURE ;
}

// tells us if we missed a record we were looking for
static int
nonexistent_record( const char *message )
{
  fprintf( stderr , "[IO] I cannot find record for tag %s \n" , message ) ;
  return FAILURE ;
}

//
int
read_propheader( struct propagator *prop )
{
  // dimensions in the file
  size_t dims[ ND ] ;
  const int MAX_HEADER_LINES = 32 ;
  char line[ MAX_LINE_LENGTH ] ;

  // error flags
  int n = 0 , dimsflag = 0 , originflag = 0 ;
  int endflag = 0 , srcflag = 0 , precflag = 0 ;
  int basisflag = 0 ;

  // loop up to some large number of possible tags
  while( n < MAX_HEADER_LINES ) {

    // read a line
    if( fgets( line , MAX_LINE_LENGTH , prop -> file ) == NULL ) {
      fprintf( stderr , "[IO] prop Header reading failed ... Leaving \n" ) ;
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
    // break when we hit the desired end_header
    if( are_equal( line , "<end_header>\n" ) ) {
      break ;
    }
    n++ ;
  }

  // check the flags have been hit
  if( dimsflag == 0 ) return nonexistent_record( "Lattice:" ) ; 
  if( originflag == 0 ) return nonexistent_record( "SrcPos:" ) ; 
  if( endflag == 0 ) return nonexistent_record( "Endian:" ) ; 
  if( precflag == 0 ) return nonexistent_record( "Precision:" ) ;
  if( srcflag == 0 ) return nonexistent_record( "Source:" ) ;
  if( basisflag == 0 ) return nonexistent_record( "Basis:" ) ;
  if( n == MAX_HEADER_LINES ) return nonexistent_record( "<end header>" ) ;

  // the NRQCD code counts from 1 instead of zero, shift to c-counting
  // instead of Fortran counting
  if( prop -> basis == NREL_FWD || prop -> basis == NREL_BWD ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      prop -> origin[ mu ] -= 1 ;
    }
  }

  // if everything is kosher we leave successfully
  return SUCCESS ;
}

// open the files and parse the header
int
read_propheaders( struct propagator *prop ,
		  const struct input_info inputs )
{
  size_t i ;
  for( i = 0 ; i < inputs.nprops ; i++ ) {
    // read and check 'em
    if( read_propheader( &prop[ i ] ) == FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
