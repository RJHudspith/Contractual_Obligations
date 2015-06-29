/**
   @file cl_diagonal.h
   @brief prototype functions for the flavour off-diagonal conserved-local contractions
 */
#ifndef CL_DIAGONAL_H
#define CL_DIAGONAL_H

/**
   @fn int cl_offdiagonal( struct propagator prop1 , struct propagator prop2 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour off-diagonal conserved-local contractions
   @return #SUCCESS or #FAILURE
 */
int
cl_offdiagonal( struct propagator prop1 , 
		struct propagator prop2 ,
		const struct site *lat ,
		const struct cut_info CUTINFO ,
		const char *outfile ) ;

#endif
