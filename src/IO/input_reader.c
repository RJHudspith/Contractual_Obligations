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

   Some of it was taken from my library GLU, some is new
*/

#include "common.h"

#include <errno.h>

#include "input_baryons.h"  // baryon contraction logic
#include "input_general.h"  // general setup information get_dims() etc..
#include "input_mesons.h"   // meson contraction logic
#include "input_tetras.h"   // tetraquark contraction logic
#include "input_VPF.h"      // VPF contraction logic
#include "input_WME.h"      // WME contraction logic

// is this what valgrind dislikes?
static struct inputs *INPUT = NULL ;

// counter for the number of tags
static size_t NTAGS = 0 ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// allocate the input file struct
static void 
pack_inputs( FILE *setup )
{
  INPUT = ( struct inputs* )calloc( MAX_TOKENS , sizeof( struct inputs ) ) ;
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
int
tag_failure( const char *tag )
{
  fprintf( stderr , "[IO] Failure looking for tag %s in input file"
	   "... Leaving\n" , tag ) ;
  return FAILURE ;
}

// look for a tag and return the index if found, else return -1
int
tag_search( const char *tag ) 
{
  size_t i ;
  for( i = 0 ; i < NTAGS ; i++ ) {
    if( are_equal( INPUT[i].TOKEN , tag ) ) return i ;
  }
  return FAILURE ;
}

//
int
unexpected_NULL( void ) 
{
  fprintf( stderr , "[IO] unexpected NULL in contraction string \n" ) ;
  return FAILURE ;
}

// get the map for our props
int
get_contraction_map( size_t *map ,
		     const char *token ,
		     const size_t nprops ) 
{
  errno = 0 ;
  char *endptr ;
  *map = (size_t)strtol( token , &endptr , 10 ) ;
  if( token == endptr ) {
    fprintf( stderr , "[IO] contraction mapping expects an integer \n" ) ;
    return FAILURE ;
  }
  if( *map >= nprops || errno == ERANGE ) {
    fprintf( stderr , "[IO] non-sensical contraction index %zu \n" , *map ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// little wrapper to free our inputs. I invisage it will grow with time
void
free_inputs( struct input_info inputs ) 
{
  free( inputs.baryons ) ;
  free( inputs.diquarks ) ;
  free( inputs.mesons ) ;
  free( inputs.tetras ) ;
  free( inputs.VPF ) ;
  free( inputs.wme ) ;
  return ;
}

// free the propagators, another structure that might get quite big
void
free_props( struct propagator *props , 
	    const size_t nprops )
{
  size_t i ;
  for( i = 0 ; i < nprops ; i++ ) {
    fclose( props[i].file ) ;
  }
  free( props ) ;
  return ;
}

// pack the input_info struct with all our contractions
int
get_input_data( struct propagator **prop ,
		struct input_info *inputs ,
		const char *file_name )
{
  // open the input file in here and free it at the bottom
  FILE *infile = fopen( file_name , "r" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[IO] input file cannot be read ... Leaving\n" ) ;
    return FAILURE ;
  }
  
  //
  int STATUS = SUCCESS ;

  // if we can open the file we push it into a big structure
  pack_inputs( infile ) ;

  // get the header
  Latt.head = header_type( INPUT ) ;
  if( Latt.head == FAILURE ) STATUS = FAILURE ;

  // read the momentum cuts struct
  if( read_cuts_struct( &( inputs -> CUTINFO ) , INPUT ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // read the dimensions from the input file
  if( get_dims( inputs -> dims , INPUT ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // initialise
  inputs -> baryons = NULL ;
  inputs -> diquarks = NULL ;
  inputs -> mesons = NULL ;
  inputs -> tetras = NULL ;
  inputs -> VPF = NULL ;
  inputs -> wme = NULL ;

  // read some props, although strange I wouldn't consider no props being an error
  get_props( *prop , &( inputs -> nprops ) , INPUT , GLU_TRUE ) ;
  *prop = (struct propagator*)malloc( ( inputs -> nprops ) * sizeof( struct propagator ) ) ;
  if( inputs -> nprops == 0 ) { 
    fprintf( stderr , "[IO] No propagator files specified \n" ) ;
    STATUS = FAILURE ;
  }
  get_props( *prop , &( inputs -> nprops ) , INPUT , GLU_FALSE ) ;

  // baryon stuff
  baryon_contractions( inputs -> baryons , &( inputs -> nbaryons ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> baryons = (struct baryon_info*)malloc( ( inputs -> nbaryons ) * sizeof( struct baryon_info ) ) ;
  if( baryon_contractions( inputs -> baryons , &( inputs -> nbaryons ) , INPUT ,
			   inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // diquark contractions
  diquark_contractions( inputs -> diquarks , &( inputs -> ndiquarks ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> diquarks = (struct meson_info*)malloc( ( inputs -> ndiquarks ) * sizeof( struct meson_info ) ) ;
  if( diquark_contractions( inputs -> diquarks , &( inputs -> ndiquarks ) , INPUT ,
			    inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // meson contractions
  meson_contractions( inputs -> mesons , &( inputs -> nmesons ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> mesons = (struct meson_info*)malloc( ( inputs -> nmesons ) * sizeof( struct meson_info ) ) ;
  if( meson_contractions( inputs -> mesons , &( inputs -> nmesons ) , INPUT ,
			  inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // tetras
  tetra_contractions( inputs -> tetras , &( inputs -> ntetras ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> tetras = (struct tetra_info*)malloc( ( inputs -> ntetras ) * sizeof( struct tetra_info ) ) ;
  if( tetra_contractions( inputs -> tetras , &( inputs -> ntetras ) , INPUT ,
			  inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // vacuum polarisation stuff
  VPF_contractions( inputs -> VPF , &( inputs -> nVPF ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> VPF = (struct VPF_info*)malloc( ( inputs -> nVPF ) * sizeof( struct VPF_info ) ) ;
  if( VPF_contractions( inputs -> VPF , &( inputs -> nVPF ) , INPUT ,
			inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // matrix element stuff
  matrix_element_contractions( inputs -> wme , &( inputs -> nWME ) , INPUT , inputs -> nprops , GLU_TRUE ) ;
  inputs -> wme = (struct WME_info*)malloc( ( inputs -> nWME ) * sizeof( struct WME_info ) ) ;
  if( matrix_element_contractions( inputs -> wme , &( inputs -> nWME ) , INPUT ,
				   inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // close the file and deallocate the buffer
  fclose( infile ) ;
  unpack_inputs( ) ;

  // if we have hit ANY problem we return FAILURE this causes it to exit
  return STATUS ;
}

