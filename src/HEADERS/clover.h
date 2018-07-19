/**
   @file clover.h
   @brief O(a^4) clover-leaf computation
 */
#ifndef CLOVER_H
#define CLOVER_H

/**
   @fn void compute_clovers( struct NRQCD_fields *F , const struct site *lat , const size_t t , const double U_0 )
   @brief compute the improved clover fields on timeslice t
 */
void
compute_clovers( struct NRQCD_fields *F ,
		 const struct site *lat ,
		 const size_t t ,
		 const double U_0 ) ;

#endif
