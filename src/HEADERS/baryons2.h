/**
   @file baryons2.h
   @brief baryon contraction code prototype functions
 */

#ifndef BARYONS2_H
#define BARYONS2_H

/**
   @fn int baryons_2diagonal( struct propagator prop , const char *outfile )
   @brief (2)flavour diagonal baryon contractions

   @return #SUCCESS or #FAILURE
 */
int
baryons_2diagonal( struct propagator prop1 ,
		   struct propagator prop2 ,
		   const char *outfile ) ;

#endif
