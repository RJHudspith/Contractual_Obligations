/**
   @file wrap_mesons.h
   @brief wrapper for the various meson contraction codes
 */

#ifndef WRAP_MESONS_H
#define WRAP_MESONS_H

/**
   @fn int contract_mesons( FILE **fprops , const struct meson_info *mesons , const int nmesons )
   @return #SUCCESS or #FAILURE
 */
int
contract_mesons( FILE **fprops ,
		 const struct meson_info *mesons ,
		 const int nmesons ) ;

#endif
