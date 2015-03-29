/**
   @fn wrap_baryons.h
   @brief prototype functions for baryon measurements
 */

#ifndef WRAP_BARYONS_H
#define WRAP_BARYONS_H

/**
   @fn int contract_baryons( struct propagator *prop , const struct meson_info *baryons , const int nbaryons )
   @brief baryon contraction logic
   @return #SUCCESS or #FAILURE
 */
int
contract_baryons( struct propagator *prop ,
		  const struct meson_info *baryons ,
		  const int nbaryons ) ;

#endif
