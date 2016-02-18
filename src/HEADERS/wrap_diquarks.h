/**
   @file wrap_diquark.h
   @brief prototype declarations for diquark wrapper
 */
#ifndef WRAP_DIQUARK_H
#define WRAP_DIQUARK_H

/**
   @fn int contract_diquarks( struct propagator *prop , const struct meson_info *diquarks , const struct cut_info CUTINFO , const size_t ndiquarks )
 */
int
contract_diquarks( struct propagator *prop ,
		   const struct meson_info *diquarks ,
		   const struct cut_info CUTINFO ,
		   const size_t ndiquarks ) ;

#endif
