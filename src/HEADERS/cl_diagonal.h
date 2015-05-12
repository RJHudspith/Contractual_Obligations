/**
   @file conserved_local.h
   @brief prototype functions for the conserved-local contractions
 */
#ifndef CONSERVED_LOCAL_H
#define CONSERVED_LOCAL_H

/**
   @fn int cl_diagonal( struct propagator prop , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour diagonal conserved-local contractions
   @return #SUCCESS or #FAILURE
*/
int
cl_diagonal( struct propagator prop ,
	     const struct site *lat ,
	     const struct cut_info CUTINFO ,
	     const char *outfile ) ;

#endif
