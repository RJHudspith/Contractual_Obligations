/**
   @file HAL_su2.h
   @brief prototype declarations for HAL driving code
 */
#ifndef HAL_SU2_H
#define HAL_SU2_H

/**
   @fn int HAL_su2( struct propagator prop1 , struct cut_info CUTINFO , const char *outfile )
   @brief performs the HAL rho-rho scattering calculation
 */
int
HAL_su2( struct propagator prop1 ,
	 struct cut_info CUTINFO ,
	 const char *outfile ) ;

#endif
