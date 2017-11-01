/**
   @file read_propheader.h
   @brief prototype functions for reading the propagator file header
 */

#ifndef READ_PROPHEADER_H
#define READ_PROPHEADER_H

/**
   @fn int read_propheader( struct propagator *prop )
   @brief read and check a propagator file header
   @return #SUCCESS or #FAILURE
 */
int
read_propheader( struct propagator *prop ) ;

/**
   @fn int read_propheaders( struct propagator *prop , const size_t nprops )
   @brief open and check our propagator files
   @param prop :: propagator file
   @return #SUCCESS or #FAILURE
 */
int
read_propheaders( struct propagator *prop ,
		  const size_t nprops ) ;

#endif
