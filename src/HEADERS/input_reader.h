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
   @fn int get_input_data( struct input_info *inputs , const char *file_name )
   @brief set the gauge field header and pass a propagator name for now
   @warning several objects in inputs get malloc'd here
 */
int
get_input_data( struct input_info *inputs , 
		const char *file_name ) ;

#endif
