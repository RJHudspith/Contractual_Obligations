/**
   @file local_local.h
   @brief prototype functions for the local-local VPF
 */

#ifndef LOCAL_LOCAL_H
#define LOCAL_LOCAL_H

/**
   @fn int local_local( FILE *prop1 , const proptype proptype1 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local flavour diagonal

   @return #SUCCESS or #FAILURE
 */
int
local_local( FILE *prop1 , 
	     const proptype proptype1 ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile ) ;

/**
   @fn int local_local_double( FILE *prop1 , const proptype proptype1 , FILE *prop2 , const proptype proptype2 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local current flavour off-diagonal

   @return #SUCCESS or #FAILURE
 */
int
local_local_double( FILE *prop1 , 
		    const proptype proptype1 ,
		    FILE *prop2 , 
		    const proptype proptype2 ,
		    const struct site *lat ,
		    const struct cut_info CUTINFO ,
		    const char *outfile ) ;

#endif
