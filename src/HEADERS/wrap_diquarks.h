/**
   @file wrap_diquarks.h
   @brief prototype declarations for diquark wrapper
 */
#ifndef WRAP_DIQUARKS_H
#define WRAP_DIQUARKS_H

/**
   @fn int contract_diquarks( struct propagator *prop , const struct meson_info *diquarks , const struct cut_info CUTINFO , const size_t ndiquarks )
   @brief driver for diquark contractions
 */
int
contract_diquarks( struct propagator *prop ,
		   const struct meson_info *diquarks ,
		   const struct cut_info CUTINFO ,
		   const size_t ndiquarks ) ;

#endif
