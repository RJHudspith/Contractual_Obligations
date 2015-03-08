/**
   @file momspace_PImunu.h
   @brief prototype functions for the momentum-space VPF
 */
#ifndef MOMSPACE_PIMUNU
#define MOMSPACE_PIMUNU

/**
   @fn void momspace_PImunu( struct PIdata *DATA_AA , struct PIdata *DATA_VV , const PImunu_projtype PROJ , const struct cut_info CUTINFO , const char *outfile )
   @brief FFTs (if available) the config-space \f$ \Pi_{\mu\nu}(q) \f$
 */
void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ) ;

#endif
