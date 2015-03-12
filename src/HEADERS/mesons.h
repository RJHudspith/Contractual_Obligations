/**
   @file mesons.h
   @brief prototype functions for meson computations
 */

#ifndef MESONS_H
#define MESONS_H

/**
   @fn int single_mesons( struct propagator prop , const char *outfile )
   @brief compute mesons from a single propagator

   @return #SUCCESS or #FAILURE
 */
int
single_mesons( struct propagator prop ,
	       const char *outfile ) ;

/**
   @fn int double_mesons( struct propagator prop1 , struct propagator prop2 , const char *outfile )
   @brief compute mesons from two propagators

   @return #SUCCESS or #FAILURE
 */
int
double_mesons( struct propagator prop1 ,
	       struct propagator prop2 ,
	       const char *outfile ) ;

#endif
