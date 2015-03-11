/**
   @file WME.h
   @brief prototype functions for Weak Matrix Elements
 */

#ifndef WEAK_MATRIX_ELEMENT_H
#define WEAK_MATRIX_ELEMENT_H

/**
   @fn int WME( FILE *S0 , const proptype s0proptype , FILE *D0 , const proptype d0proptype , FILE *S1 , const proptype s1proptype , FILE *D1 , const proptype d1proptype , const char *outfile )
   @brief weak matrix element four quark operator insertion
   @param S0 :: strange quark propagator at Wall 0
   @param D0 :: down quark propagator at Wall 0
   @param S1 :: strange quark propagator at Wall LT/2
   @param D1 :: down quark propagator at Wall LT/2

   @warning propagators have to be wall sources

   @return #FAILURE or #SUCCESS
 */
int
WME( FILE *S0 , const proptype s0proptype ,
     FILE *D0 , const proptype d0proptype ,
     FILE *S1 , const proptype s1proptype ,
     FILE *D1 , const proptype d1proptype ,
     const char *outfile ) ;

#endif
