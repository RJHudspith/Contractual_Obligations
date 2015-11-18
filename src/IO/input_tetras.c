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
	      const int nprops ,
	      const int meas_idx ,
	      const char *message ) 
{
  printf( "\n" ) ;
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
  printf( "[IO] %s_%d :: Contracting prop %d with prop %d with prop %d\n" , 
	  message , meas_idx , tetras -> map[0] , tetras -> map[1] , tetras -> map[2] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( tetras -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%d :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%d :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// tetra expects four contraction indices
int
tetra_contractions( struct tetra_info *tetras , 
		    int *ntetras ,
		    struct inputs *INPUT ,
		    const int nprops ,
		    const GLU_bool first_pass ) 
{
  *ntetras = 0 ;
  char str[ 32 ] ;
  while( *ntetras < MAX_CONTRACTIONS ) {
    sprintf( str , "TETRA%d" , *ntetras ) ;
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
