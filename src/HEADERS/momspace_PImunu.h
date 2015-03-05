/**
   @file momspace_PImunu.h
   @brief prototype functions for the momentum-space VPF
 */
#ifndef MOMSPACE_PIMUNU
#define MOMSPACE_PIMUNU

/**
   @fn void momspace_PImunu( struct PIdata *DATA_AA , struct PIdata *DATA_VV )
   @brief FFTs (if available) the config-space \f$ \Pi_{\mu\nu}(q) \f$
 */
void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV ) ;

#endif
