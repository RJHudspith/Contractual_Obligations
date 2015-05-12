/**
   @file ll_diagonal.h
   @brief prototype functions for the local-local VPF
 */
#ifndef LL_DIAGONAL_H
#define LL_DIAGONAL_H

/**
   @fn int ll_diagonal( struct propagator prop , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local flavour diagonal

   @return #SUCCESS or #FAILURE
 */
int
ll_diagonal( struct propagator prop ,
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
