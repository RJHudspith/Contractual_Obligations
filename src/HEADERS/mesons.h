/**
   @file mesons.h
   @brief prototype functions for meson computations
 */

#ifndef MESONS_H
#define MESONS_H

/**
   @fn int single_mesons( FILE *fprop , const proptype proptype1 , const char *outfile )
   @brief compute mesons from a single propagator

   @return #SUCCESS or #FAILURE
 */
int
single_mesons( FILE *fprop ,
	       const proptype proptype1 ,
	       const char *outfile ) ;

/**
   @fn int double_mesons( FILE *prop1 , const proptype proptype1 , FILE *prop2 , const proptype proptype2 , const char *outfile )
   @brief compute mesons from two propagators

   @return #SUCCESS or #FAILURE
 */
int
double_mesons( FILE *prop1 , 
	       const proptype proptype1 ,
	       FILE *prop2 ,
	       const proptype proptype2 ,
	       const char *outfile ) ;

#endif
