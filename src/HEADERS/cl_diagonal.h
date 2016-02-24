/**
   @file cl_diagonal.h
   @brief prototype functions for the conserved-local contractions
 */
#ifndef CL_DIAGONAL_H
#define CL_DIAGONAL_H

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
