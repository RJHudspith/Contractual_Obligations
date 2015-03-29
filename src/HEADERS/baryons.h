/**
   @file baryons.h
   @brief baryon contraction code prototype functions
 */

#ifndef BARYONS_H
#define BARYONS_H

/**
   @fn int baryons_diagonal( struct propagator prop , const char *outfile )
   @brief flavour diagonal baryon contractions
   @return #SUCCESS or #FAILURE
 */
int
baryons_diagonal( struct propagator prop ,
		  const char *outfile ) ;

#endif
