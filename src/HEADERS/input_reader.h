/**
   @file input_reader.h
   @brief read the input file

   Expects an input file of the form

   TAG = VALUE
 */

#ifndef INPUT_READER_H
#define INPUT_READER_H

/**
   @fn void free_inputs( struct input_info inputs )
   @brief free the input struct
 */
void
free_inputs( struct input_info inputs ) ;

/**
   @fn void free_props( struct propagator *prop , const size_t nprops )
   @brief free the propagator storage
 */
void
free_props( struct propagator *prop , 
	    const size_t nprops ) ;

/**
   @fn int get_contraction_map( int *map , const char *token , const size_t nprops ) 
   @brief contraction map for our propagators
 */
int
get_contraction_map( size_t *map ,
		     const char *token ,
		     const size_t nprops ) ;

/**
   @fn int get_input_data( struct propagator **prop , struct input_info *inputs , const char *file_name )
   @brief set the gauge field header and pass a propagator name for now
   @warning several objects in inputs get malloc'd here
 */
int
get_input_data( struct propagator **prop ,  
		struct input_info *inputs , 
		const char *file_name ) ;

/**
   @fn int tag_failure( const char *tag )
   @brief prints out a problem tag
 */
int
tag_failure( const char *tag ) ;

/**
   @fn int tag_search( const char *tag )
   @brief search for a tag
 */
int
tag_search( const char *tag ) ;

/**
   @fn int unexpected_NULL( void )
   @brief error message for an unexpected NULL
 */
int
unexpected_NULL( void ) ;

#endif
