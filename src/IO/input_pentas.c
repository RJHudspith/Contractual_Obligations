/**
   @file input_pentas.c
   @brief pentaquark contractions dictated by the input file
 */
#include "common.h"

#include "input_reader.h"

// tokenize the map of pentas
static int
penta_tokens( struct penta_info *pentas ,
	      const char *penta_str ,
	      const size_t nprops ,
	      const size_t meas_idx ,
	      const char *message )
{
  fprintf( stdout , "\n" ) ;
  // starts with the contraction indices
  char *token = (char*)strtok( (char*)penta_str , "," ) ;
  if( get_contraction_map( &( pentas -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( pentas -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( pentas -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  fprintf( stdout , "[IO] %s_%zu :: Contracting prop %zu with prop %zu "
	   "with prop %zu \n" , 
	   message , meas_idx , pentas -> map[0] , pentas -> map[1] , 
	   pentas -> map[2]  ) ;
    // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( pentas -> outfile , "%s" , token ) ;
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

// penta expects three contraction indices
int
penta_contractions( struct penta_info *pentas , 
		    size_t *npentas ,
		    struct inputs *INPUT ,
		    const size_t nprops ,
		    const GLU_bool first_pass ) 
{
  *npentas = 0 ;
  char str[ 32 ] ;
  while( *npentas < MAX_CONTRACTIONS ) {
    sprintf( str , "PENTA%zu" , *npentas ) ;
    const int penta_idx = tag_search( str ) ;
    if( penta_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( penta_tokens( &pentas[ *npentas ] , 
			INPUT[ penta_idx ].VALUE , 
			nprops , *npentas , "Penta" ) 
	  == FAILURE ) return FAILURE ;
    }
    *npentas = *npentas + 1 ;
  }
  return SUCCESS ;
}
