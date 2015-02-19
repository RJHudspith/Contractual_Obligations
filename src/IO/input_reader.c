/*
    Copyright 2013 Renwick James Hudspith

    This file (input_reader.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file input_reader.c
   @brief reader for the input file
*/

#include "common.h"

// we might want to change this at some point
#define GLU_STR_LENGTH 96

// maximum number of tokens
#define INPUTS_LENGTH 36

// tokenize the input file
struct inputs {
  char TOKEN[ GLU_STR_LENGTH ] ;
  char VALUE[ GLU_STR_LENGTH ] ;
} ;

// is this what valgrind dislikes?
static struct inputs *INPUT ;

// counter for the number of tags
static int NTAGS = 0 ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// allocate the input file struct
static void 
pack_inputs( FILE *setup )
{
  INPUT = ( struct inputs* )malloc( INPUTS_LENGTH * sizeof( struct inputs ) ) ;
  // and put into the buffer
  while( NTAGS++ , fscanf( setup , "%s = %s" , INPUT[ NTAGS ].TOKEN , INPUT[ NTAGS ].VALUE )  != EOF ) { }
  return ;
}

// frees up the struct
static void
unpack_inputs( void )
{
  free( INPUT ) ;
  return ;
}

// prints out the problematic tag
static int
tag_failure( const char *tag )
{
  printf( "[IO] Failure looking for tag %s in input file ... Leaving\n" , tag ) ;
  return FAILURE ;
}

// look for a tag and return the index if found, else return -1
static int
tag_search( const char *tag ) 
{
  int i ;
  for( i = 0 ; i < NTAGS ; i++ ) {
    if( are_equal( INPUT[i].TOKEN , tag ) ) return i ;
  }
  return FAILURE ;
}

// quickly get the configuration number from the input file
static int
confno( void )
{
  return atoi( INPUT[tag_search( "CONFNO" )].VALUE ) ;
}

// get the header type
static header_mode
header_type( void )
{
  const int header_idx = tag_search( "HEADER" ) ;
  if( header_idx == FAILURE ) { return tag_failure( "HEADER" ) ; }

  if( are_equal( INPUT[header_idx].VALUE , "NERSC" ) ) {
    printf( "[IO] Attempting to read a NERSC file \n" ) ;
    return NERSC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "HIREP" ) ) {
    printf( "[IO] Attempting to read a HIREP file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return HIREP_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "MILC" ) ) {
    printf( "[IO] Attempting to read a MILC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return MILC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "SCIDAC" ) ) {
    printf( "[IO] Attempting to read a SCIDAC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "LIME" ) ) {
    printf( "[IO] Attempting to read an LIME file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    printf( "[IO] WARNING!! NOT CHECKING ANY CHECKSUMS!! \n" ) ;
    return LIME_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_SCIDAC" ) ) {
    printf( "[IO] Attempting to read an ILDG (Scidac) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return ILDG_SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_BQCD" ) ) {
    printf( "[IO] Attempting to read an ILDG (BQCD) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return ILDG_BQCD_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "UNIT" ) ) {
    printf( "[IO] Generating a UNIT gauge configuration \n" ) ;
    return UNIT_GAUGE ;
  }
  printf( "[IO] HEADER %s not recognised ... Leaving \n" , 
	  INPUT[header_idx].VALUE ) ;
  return FAILURE ; 
}

// fills the INFILE struct with all the useful information
int
get_input_data( char prop[][ GLU_STR_LENGTH ] ,
		int *nprops , 
		int *dims ,
		const char *file_name )
{
  // open the input file in here and free it at the bottom
  FILE *infile = fopen( file_name , "r" ) ;
  if( infile == NULL ) {
    printf( "[IO] input file cannot be read ... Leaving\n" ) ;
    return FAILURE ;
  }
  
  //
  int STATUS = SUCCESS ;

  // if we can open the file we push it into a big structure
  pack_inputs( infile ) ;

  int INPUT_FAILS = 0 ; // counter for the number of failures

  // get the header
  Latt.head = header_type( ) ;
  if( Latt.head == FAILURE ) STATUS = FAILURE ;

  // read some props
  int count = 0 ;
  char str[ 32 ] ;
  while( count < *nprops ) {
    sprintf( str , "PROP%d" , count ) ;
    const int prop_idx = tag_search( str ) ;
    if( prop_idx == FAILURE ) break ;
    strcpy( prop[ count ] , INPUT[ prop_idx ].VALUE ) ;
    count++ ;
  }
  if( count == 0 ) { 
    printf( "[IO] No propagator files specified ... Leaving \n" ) ;
    STATUS = FAILURE ;
  }
  *nprops = count ;

  // read the dimensions from the input file
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    sprintf( str , "DIMS%d" , mu ) ;
    const int dims_idx = tag_search( str ) ;
    if( dims_idx == FAILURE ) { 
      tag_failure( str ) ; 
      STATUS = FAILURE ; 
      break ; 
    }
    dims[ mu ] = (int)atoi( INPUT[ dims_idx ].VALUE ) ;
  }

  // close the file and deallocate the buffer
  fclose( infile ) ;
  unpack_inputs( ) ;

  // if we have hit ANY problem we return GLU_FAILURE this causes it to exit
  return STATUS ;
}
