/**
   @file wrap_mesons.h
   @brief wrapper for the various meson contraction codes
 */

#ifndef WRAP_MESONS_H
#define WRAP_MESONS_H

/**
   @fn int contract_mesons( struct propagator *prop , const struct meson_info *mesons , const struct cut_info CUTINFO , const size_t nmesons )
   @return #SUCCESS or #FAILURE
 */
int
contract_mesons( struct propagator *prop ,
		 const struct meson_info *mesons ,
		 const struct cut_info CUTINFO ,
		 const size_t nmesons ) ;

#endif
