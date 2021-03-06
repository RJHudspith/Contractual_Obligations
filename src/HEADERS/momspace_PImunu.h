/**
   @file momspace_PImunu.h
   @brief prototype functions for the momentum-space VPF
 */
#ifndef MOMSPACE_PIMUNU
#define MOMSPACE_PIMUNU

/**
   @fn void momspace_PImunu( struct PIdata *AA , struct PIdata *VV , const struct cut_info CUTINFO , const char *outfile , const current_type current )
   @brief FFTs (if available) the config-space \f$ \Pi_{\mu\nu}(q) \f$

   @warning overwrites the PIdatas with their fourier transform
 */
void
momspace_PImunu( struct PIdata *AA ,
		 struct PIdata *VV ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ,
		 const current_type current ) ;

#endif
