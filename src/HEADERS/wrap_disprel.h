/**
   @file wrap_disprel.h
   @brief prototype declarations for dispersion relation wrapper
 */
#ifndef WRAP_DISPREL_H
#define WRAP_DISPREL_H

/**
   @fn int contract_disprels( struct propagator *prop , const struct dispersion_info *dispersions , const struct cut_info CUTINFO , const int ndispersions )
   @brief dispersion relation contraction wrapper
   @return #SUCCESS or #FAILURE
 */
int
contract_disprels( struct propagator *prop ,
		   const struct dispersion_info *dispersions ,
		   const struct cut_info CUTINFO ,
		   const int ndispersions ) ;

#endif
