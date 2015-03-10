/**
   @file read_propheader.h
   @brief prototype functions for reading the propagator file header
 */

#ifndef READ_PROPHEADER_H
#define READ_PROPHEADER_H

/**
   @fn int read_check_header( FILE *propfile , const GLU_bool first_read )
   @brief read and check a propagator file header
   @return #SUCCESS or #FAILURE
 */
int
read_check_header( FILE *propfile ,
		   const GLU_bool first_read ) ;

/**
   @file int read_propheaders( FILE **fprops , const struct input_info inputs )
   @brief open and check our propagator files
   @return #SUCCESS or #FAILURE
 */
int
read_propheaders( FILE **fprops ,
		  const struct input_info inputs ) ;

#endif
