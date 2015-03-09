/**
   @file input_reader.h
   @brief read the input file

   Expects an input file of the form

   TAG = VALUE
 */

#ifndef INPUT_READER_H
#define INPUT_READER_H

/**
   @fn int
   @brief set the gauge field header and pass a propagator name for now
 */
int
get_input_data( char prop[][GLU_STR_LENGTH] ,
		int *nprops ,
		struct meson_info *mesons ,
		int *nmesons ,
		struct VPF_info *VPF ,
		int *nVPF ,
		struct cut_info *CUTINFO ,
		int *dims ,
		const char *file_name ) ;

#endif
