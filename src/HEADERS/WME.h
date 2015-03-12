/**
   @file WME.h
   @brief prototype functions for Weak Matrix Elements
 */

#ifndef WEAK_MATRIX_ELEMENT_H
#define WEAK_MATRIX_ELEMENT_H

/**
   @fn int WME( struct propagator s0 , struct propagator d0 , struct propagator s1 , struct propagator d1 , const char *outfile )
   @brief weak matrix element four quark operator insertion
   @param s0 :: strange quark propagator at Wall 0
   @param d0 :: down quark propagator at Wall 0
   @param s1 :: strange quark propagator at Wall LT/2
   @param d1 :: down quark propagator at Wall LT/2
   @param outfile :: output file
   @warning propagators have to be wall sources

   @return #FAILURE or #SUCCESS
 */
int
WME( struct propagator s0 ,
     struct propagator d0 ,
     struct propagator s1 ,
     struct propagator d1 ,
     const char *outfile ) ;

#endif
