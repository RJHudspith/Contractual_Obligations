/**
   @file ll_offdiagonal.h
   @brief prototype functions for the flavour off diagonal local-local VPF
 */
#ifndef LL_OFFDIAGONAL_H
#define LL_OFFDIAGONAL_H

/**
   @fn int ll_offdiagonal( struct propagator prop1 , struct propagator prop2 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief local-local current flavour off-diagonal
   @return #SUCCESS or #FAILURE
 */
int
ll_offdiagonal( struct propagator prop1 ,
		struct propagator prop2 ,
		const struct site *lat ,
		const struct cut_info CUTINFO ,
		const char *outfile ) ;

#endif
