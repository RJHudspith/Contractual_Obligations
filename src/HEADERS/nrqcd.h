/**
   @file nrqcd.h
   @brief on the fly NRQCD code prototype declarations
 */
#ifndef NRQCD_H
#define NRQCD_H

/**
   @fn int compute_nrqcd_props( struct propagator *prop , const size_t nprops ) 
   @brief computes all of the NRQCD propagators we ask for
 */
int
compute_nrqcd_props( struct propagator *prop ,
		     const size_t nprops ) ;

/**
   @fn int free_nrqcd_props( struct propagator *prop , const size_t nprops )
   @brief free the allocated NRQCD propagators
 */
int
free_nrqcd_props( struct propagator *prop ,
		  const size_t nprops ) ;

#endif
