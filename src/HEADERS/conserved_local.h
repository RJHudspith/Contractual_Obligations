/**
   @file conserved_local.h
   @brief prototype functions for the conserved-local contractions
 */
#ifndef CONSERVED_LOCAL_H
#define CONSERVED_LOCAL_H

/**
   @fn int conserved_local( FILE *prop1 , const proptype proptype1 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour diagonal conserved-local contractions
   @return #SUCCESS or #FAILURE
*/
int
conserved_local( FILE *prop1 , 
		 const proptype proptype1 , 
		 const struct site *lat ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ) ;

/**
   @fn int conserved_local_double( FILE *prop1 , const proptype proptype1 , FILE *prop2 , const proptype proptype2 , const struct site *lat , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour off-diagonal conserved-local contractions
   @return #SUCCESS or #FAILURE
 */
int
conserved_local_double( FILE *prop1 , 
			const proptype proptype1 ,
			FILE *prop2 , 
			const proptype proptype2 ,
			const struct site *lat ,
			const struct cut_info CUTINFO ,
			const char *outfile ) ;
#endif
