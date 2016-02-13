/**
   @file input_VPF.c
   @brief VPF contraction logic from the input file
 */
#include "common.h"

#include "input_reader.h"

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// get the various contraction types
static int
get_current_type( current_type *current ,
		  const char *token ,
		  const size_t meas ) 
{
  if( are_equal( token , "LOCAL_LOCAL" ) ) {
    *current = LOCAL_LOCAL ;
  } else if( are_equal( token , "CONSERVED_LOCAL" ) ) {
    *current = CONSERVED_LOCAL ;
  } else {
    printf( "[IO] I don't understand VPF type %s\n" , token ) ;
    return FAILURE ;
  }
  printf( "[IO] VPF_%zu :: We are performing %s contractions \n" , 
	  meas , token ) ;
  return SUCCESS ;
}

// tokenize each vpf contraction
static int
VPF_tokens( struct VPF_info *VPF ,
	    const char *VPF_str ,
	    const size_t nprops ,
	    const size_t meas_idx ) 
{
  printf( "\n" ) ;

  // starts with the propagator map
  char *token = (char*)strtok( (char*)VPF_str , "," ) ;
  if( get_contraction_map( &( VPF -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( VPF -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] VPF_%zu :: Contracting prop %zu with prop %zu \n" , 
	  meas_idx , VPF -> map[0] , VPF -> map[1] ) ;

  // check for current
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_current_type( &( VPF -> current ) , token , meas_idx ) == FAILURE ) return FAILURE ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( VPF -> outfile , "%s" , token ) ;
  printf( "[IO] VPF_%zu :: Contraction file in %s \n" , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] VPF_%zu :: Unexpected extra contraction info %s \n" ,
	    meas_idx , token ) ;
  }
  return SUCCESS ;
}

// meson contractions of the form prop_idx1,prop_idx2,proptype,outfile
int
VPF_contractions( struct VPF_info *VPF , 
		  size_t *nVPF ,
		  const struct inputs *INPUT ,
		  const size_t nprops ,
		  const GLU_bool first_pass ) 
{
  *nVPF = 0 ;
  char str[ 32 ] ;
  while( *nVPF < MAX_CONTRACTIONS ) {
    sprintf( str , "VPF%zu" , *nVPF ) ;
    const int VPF_idx = tag_search( str ) ;
    if( VPF_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( VPF_tokens( &VPF[ *nVPF ] , 
		      INPUT[ VPF_idx ].VALUE , 
		      nprops , *nVPF ) 
	  == FAILURE ) return FAILURE ;
    }
    *nVPF = *nVPF + 1 ;
  }
  return SUCCESS ;
}
