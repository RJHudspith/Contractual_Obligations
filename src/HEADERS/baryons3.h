/**
   @file baryons3.h
   @brief baryon contraction code prototype functions
 */

#ifndef BARYONS3_H
#define BARYONS3_H

/**
   @fn int baryons_2diagonal( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , const char *outfile )
   @brief (3)-flavour diagonal baryon contractions

   @return #SUCCESS or #FAILURE
 */
int
baryons_3fdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    struct propagator prop3 ,
		    const char *outfile ) ;

#endif
