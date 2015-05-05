/**
   @file input_disprel.c
   @brief dispersion relation contraction logic
 */
#include "common.h"

#include "input_reader.h"

// tokenize each meson-dispersion contraction, will need to think about these some more ...
static int
dispersion_tokens( struct dispersion_info *dispersions ,
		   const char *dispersion_str ,
		   const int nprops ,
		   const int meas_idx ,
		   const char *message ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)dispersion_str , "," ) ;
  if( get_contraction_map( &( dispersions -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( dispersions -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] %s_%d :: Contracting prop %d with prop %d \n" , 
	  message , meas_idx , dispersions -> map[0] , dispersions -> map[1] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( dispersions -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%d :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%d :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// dispersion contractions of the form prop_idx1,prop_idx2,proptype,outfile
int
dispersion_contractions( struct dispersion_info *dispersions , 
			 int *ndispersions ,
			 const struct inputs *INPUT ,
			 const int nprops ,
			 const GLU_bool first_pass ) 
{
  *ndispersions = 0 ;
  char str[ 32 ] ;
  while( *ndispersions < MAX_CONTRACTIONS ) {
    sprintf( str , "DISPERSION%d" , *ndispersions ) ;
    const int dispersion_idx = tag_search( str ) ;
    if( dispersion_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( dispersion_tokens( &dispersions[ *ndispersions ] , 
			     INPUT[ dispersion_idx ].VALUE , 
			     nprops , *ndispersions , "Dispersion" ) 
	  == FAILURE ) return FAILURE ;
    }
    *ndispersions = *ndispersions + 1 ;
  }
  return SUCCESS ;
}
