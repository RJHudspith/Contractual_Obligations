/**
   @file wrap_WME.h
   @brief wrapper for WME contractions
 */

#ifndef WRAP_WME_H
#define WRAP_WME_H

/**
   @fn int contract_WME( FILE **fprops , const struct WME_info *wme , const int nWME )
   @brief performs the contractions specified in the input file

   @warning propagator files must be wall sources

   @return #SUCCESS or #FAILURE
 */
int
contract_WME( FILE **fprops ,
	      const struct WME_info *wme ,
	      const int nWME ) ;

#endif
