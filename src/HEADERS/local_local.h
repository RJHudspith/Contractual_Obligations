/**
   @file local_local.h
   @brief prototype functions for the local-local VPF
 */

#ifndef LOCAL_LOCAL_H
#define LOCAL_LOCAL_H

/**
   @fn int local_local( struct propagator prop , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local flavour diagonal

   @return #SUCCESS or #FAILURE
 */
int
local_local( struct propagator prop ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile ) ;

/**
   @fn int local_local_double( struct propagator prop1 , struct propagator prop2 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local current flavour off-diagonal

   @return #SUCCESS or #FAILURE
 */
int
local_local_double( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const struct site *lat ,
		    const struct cut_info CUTINFO ,
		    const char *outfile ) ;

#endif
