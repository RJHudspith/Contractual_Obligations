/**
   @file input_mesons.c
   @brief meson contraction input file logic
 */
#include "common.h"

#include "input_reader.h"

// tokenize each meson
static int
twopoint_tokens( struct meson_info *mesons ,
		 const char *meson_str ,
		 const size_t nprops ,
		 const size_t meas_idx ,
		 const char *message ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)meson_str , "," ) ;
  if( get_contraction_map( &( mesons -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( mesons -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] %s_%zu :: Contracting prop %zu with prop %zu \n" , 
	  message , meas_idx , mesons -> map[0] , mesons -> map[1] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( mesons -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%zu :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%zu :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// diquark contractions of the form prop_idx1,prop_idx2,outfile
int
diquark_contractions( struct meson_info *diquarks , 
		      size_t *ndiquarks ,
		      const struct inputs *INPUT ,
		      const size_t nprops ,
		      const GLU_bool first_pass ) 
{
  *ndiquarks = 0 ;
  char str[ 32 ] ;
  while( *ndiquarks < MAX_CONTRACTIONS ) {
    sprintf( str , "DIQUARK%zu" , *ndiquarks ) ;
    const int diquark_idx = tag_search( str ) ;
    if( diquark_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( twopoint_tokens( &diquarks[ *ndiquarks ] , 
			   INPUT[ diquark_idx ].VALUE , 
			   nprops , *ndiquarks , "Diquark" ) 
	  == FAILURE ) return FAILURE ;
    }
    *ndiquarks = *ndiquarks + 1 ;
  }
  return SUCCESS ;
}

// meson contractions of the form prop_idx1,prop_idx2,outfile
int
meson_contractions( struct meson_info *mesons , 
		    size_t *nmesons ,
		    const struct inputs *INPUT ,
		    const size_t nprops ,
		    const GLU_bool first_pass ) 
{
  *nmesons = 0 ;
  char str[ 32 ] ;
  while( *nmesons < MAX_CONTRACTIONS ) {
    sprintf( str , "MESON%zu" , *nmesons ) ;
    const int meson_idx = tag_search( str ) ;
    if( meson_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( twopoint_tokens( &mesons[ *nmesons ] , 
			   INPUT[ meson_idx ].VALUE , 
			   nprops , *nmesons , "Meson" ) 
	  == FAILURE ) return FAILURE ;
    }
    *nmesons = *nmesons + 1 ;
  }
  return SUCCESS ;
}
