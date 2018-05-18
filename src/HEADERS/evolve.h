/**
   @file evolve.h
   @brief prototype declarations for the NRQCD evolution
 */
#ifndef EVOLVE_H
#define EVOLVE_H

/**
   @fn void compute_props( struct propagator *prop , struct NRQCD_fields *F , const struct site *lat , const size_t nprops , const double tadref )
   @brief computes all the NRQCD propagators we will use
 */
void
compute_props( struct propagator *prop ,
	       struct NRQCD_fields *F ,
	       const struct site *lat ,
	       const size_t nprops ,
	       const double tadref ) ;

#endif
