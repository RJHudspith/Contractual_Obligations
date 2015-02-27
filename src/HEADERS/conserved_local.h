/**
   @file conserved_local.h
   @brief prototype functions for the conserved-local contractions
 */
#ifndef CONSERVED_LOCAL_H
#define CONSERVED_LOCAL_H

/**
   @fn int conserved_local( FILE *fprop1 , const proptype proptype1 , const struct site *lat )
   @brief conserved-local contractions
   @return #SUCCESS or #FAILURE
*/
int
conserved_local( FILE *fprop1 , 
		 const proptype proptype1 , 
		 const struct site *lat ) ;

#endif
