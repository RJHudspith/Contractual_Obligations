/**
   @file input_WME.c
   @brief WME contraction logic from the input file
 */
#include "common.h"

#include "input_reader.h"

// tokenize each WME
static int
matrix_element_tokens( struct WME_info *wme ,
		       const char *WME_str ,
		       const int nprops ,
		       const int meas_idx ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)WME_str , "," ) ;
  if( get_contraction_map( &( wme -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[3] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] WME_%d :: Contracting prop %d with prop %d with %d with %d \n" , 
	  meas_idx , wme -> map[0] , wme -> map[1] , wme -> map[2] , wme -> map[3] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( wme -> outfile , "%s" , token ) ;
  printf( "[IO] WME_%d :: Contraction file in %s \n" , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] WME_%d :: Unexpected extra contraction info %s \n" ,
	    meas_idx , token ) ;
  }
  return SUCCESS ;
}

// matrix elemen contractions expect prop map of 4 indices
int
matrix_element_contractions( struct WME_info *wme , 
			     int *nWME ,
			     const struct inputs *INPUT ,
			     const int nprops ,
			     const GLU_bool first_pass ) 
{
  *nWME = 0 ;
  char str[ 32 ] ;
  while( *nWME < MAX_CONTRACTIONS ) {
    sprintf( str , "WME%d" , *nWME ) ;
    const int WME_idx = tag_search( str ) ;
    if( WME_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( matrix_element_tokens( &wme[ *nWME ] , 
				 INPUT[ WME_idx ].VALUE , 
				 nprops , *nWME ) 
	  == FAILURE ) return FAILURE ;
    }
    *nWME = *nWME + 1 ;
  }
  return SUCCESS ;
}
