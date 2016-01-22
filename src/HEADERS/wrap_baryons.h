/**
   @file wrap_baryons.h
   @brief prototype functions for baryon measurements
 */
#ifndef WRAP_BARYONS_H
#define WRAP_BARYONS_H

/**
   @fn int contract_baryons( struct propagator *prop , const struct baryon_info *baryons , const struct cut_info CUTINFO , const size_t nbaryons )
   @brief baryon contraction logic
   @return #SUCCESS or #FAILURE
 */
int
contract_baryons( struct propagator *prop ,
		  const struct baryon_info *baryons ,
		  const struct cut_info CUTINFO , 
		  const size_t nbaryons ) ;

#endif
