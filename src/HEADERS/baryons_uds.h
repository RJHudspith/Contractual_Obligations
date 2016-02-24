/**
   @file baryons_uds.h
   @brief (uds) baryon contraction code prototype functions
 */
#ifndef BARYONS_UDS_H
#define BARYONS_UDS_H

/**
   @fn int baryons_3fdiagonal( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , const struct cut_info CUTINFO , const char *outfile )
   @brief (3)-flavour diagonal baryon contractions

   @return #SUCCESS or #FAILURE
 */
int
baryons_3fdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    struct propagator prop3 ,
		    const struct cut_info CUTINFO ,
		    const char *outfile ) ;

#endif
