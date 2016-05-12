/**
   @file input_tetras.c
   @brief tetraquark contractions dictated by the input file
 */
#include "common.h"

#include "input_reader.h"

// tokenize the map of tetras
static int
tetra_tokens( struct tetra_info *tetras ,
	      const char *tetra_str ,
	      const size_t nprops ,
	      const size_t meas_idx ,
	      const char *message ) 
{
  fprintf( stdout , "\n" ) ;
  // starts with the contraction indices
  char *token = (char*)strtok( (char*)tetra_str , "," ) ;
  if( get_contraction_map( &( tetras -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( tetras -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( tetras -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( tetras -> map[3] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  fprintf( stdout , "[IO] %s_%zu :: Contracting prop %zu with prop %zu "
	   "with prop %zu with prop %zu\n" , 
	   message , meas_idx , tetras -> map[0] , tetras -> map[1] , 
	   tetras -> map[2] , tetras -> map[3] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( tetras -> outfile , "%s" , token ) ;
  fprintf( stdout , "[IO] %s_%zu :: Contraction file in %s \n" , 
	   message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    fprintf( stderr , "[IO] %s_%zu :: Unexpected extra contraction info %s \n" ,
	     message , meas_idx , token ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// tetra expects four contraction indices
int
tetra_contractions( struct tetra_info *tetras , 
		    size_t *ntetras ,
		    struct inputs *INPUT ,
		    const size_t nprops ,
		    const GLU_bool first_pass ) 
{
  *ntetras = 0 ;
  char str[ 32 ] ;
  while( *ntetras < MAX_CONTRACTIONS ) {
    sprintf( str , "TETRA%zu" , *ntetras ) ;
    const int tetra_idx = tag_search( str ) ;
    if( tetra_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( tetra_tokens( &tetras[ *ntetras ] , 
			INPUT[ tetra_idx ].VALUE , 
			nprops , *ntetras , "Tetra" ) 
	  == FAILURE ) return FAILURE ;
    }
    *ntetras = *ntetras + 1 ;
  }
  return SUCCESS ;
}
