/**
   @file input_baryons.c
   @brief baryon contraction logic
 */
#include "common.h"

#include "input_reader.h"

// tokenize the map of baryons
static int
threepoint_tokens( struct baryon_info *baryons ,
		   const char *baryon_str ,
		   const size_t nprops ,
		   const size_t meas_idx ,
		   const char *message ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)baryon_str , "," ) ;
  if( get_contraction_map( &( baryons -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( baryons -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( baryons -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] %s_%zu :: Contracting prop %zu with prop %zu with prop %zu \n" , 
	  message , meas_idx , baryons -> map[0] , 
	  baryons -> map[1] , baryons -> map[2] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( baryons -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%zu :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%zu :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// baryon contractions of the form prop_idx1,prop_idx2,proptype,outfile
int
baryon_contractions( struct baryon_info *baryons , 
		     size_t *nbaryons ,
		     const struct inputs *INPUT ,
		     const size_t nprops ,
		     const GLU_bool first_pass ) 
{
  *nbaryons = 0 ;
  char str[ 32 ] ;
  while( *nbaryons < MAX_CONTRACTIONS ) {
    sprintf( str , "BARYON%zu" , *nbaryons ) ;
    const int baryon_idx = tag_search( str ) ;
    if( baryon_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( threepoint_tokens( &baryons[ *nbaryons ] , 
			     INPUT[ baryon_idx ].VALUE , 
			     nprops , *nbaryons , "Baryon" ) 
	  == FAILURE ) return FAILURE ;
    }
    *nbaryons = *nbaryons + 1 ;
  }
  return SUCCESS ;
}
