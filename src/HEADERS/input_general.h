/**
   @file input_general.h
   @brief prototype declarations for general input file stuff
 */
#ifndef INPUT_GENERAL_H
#define INPUT_GENERAL_H

/**
   @fn int confno( const struct inputs *INPUT )
   @brief read the configuration number from the input file
   @return the configuration number defined in the input file
 */
int
confno( const struct inputs *INPUT ) ;

/**
   @fn int get_dims( size_t *dims , const struct inputs *INPUT )
   @brief read the lattice dimensions from the input file to compare with the propagator read
   @return #SUCCESS or #FAILURE
 */
int
get_dims( size_t *dims , 
	  const struct inputs *INPUT ) ;

/**
   @fn int get_props( struct propagator *props , size_t *nprops , const struct inputs *INPUT , const GLU_bool first_pass )
   @brief get the list of propagators we will be using in this run
   @return #SUCCESS or #FAILURE
 */
int
get_props( struct propagator *props ,
	   size_t *nprops ,
	   const struct inputs *INPUT ,
	   const GLU_bool first_pass ) ;

/**
   @fn header_mode header_type( const struct inputs *INPUT )
   @brief get the gauge configuration header type from the input file
 */
header_mode
header_type( const struct inputs *INPUT ) ;

/**
   @fn int read_cuts_struct( struct cut_info *CUTINFO , const struct inputs *INPUT )
   @brief pack the cut information struct
   @return #SUCCESS or #FAILURE
 */
int
read_cuts_struct( struct cut_info *CUTINFO ,
		  const struct inputs *INPUT ) ;

#endif
