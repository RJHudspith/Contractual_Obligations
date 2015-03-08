/**
   @file input_reader.h
   @brief read the input file

   Expects an input file of the form

   TAG = VALUE
 */

#ifndef INPUT_READER_H
#define INPUT_READER_H

/**
   @fn int get_input_data( char **prop , struct meson_info *mesons , struct cut_info *CUTINFO , int *nmesons , int *nprops , const char *file_name )
   @brief set the gauge field header and pass a propagator name for now
 */
int
get_input_data( char prop[][GLU_STR_LENGTH] ,
		struct meson_info *mesons ,
		struct cut_info *CUTINFO ,
		int *nmesons ,
		int *nprops ,
		int *dims ,
		const char *file_name ) ;

#endif
