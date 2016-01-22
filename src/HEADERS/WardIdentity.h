/**
   @file WardIdentity.h
   @brief prototype functions for ward identity checking
 */

#ifndef WARDIDENTITY_H
#define WARDIDENTITY_H

/**
   @fn void compute_p_psq( double **p , double *psq , const struct veclist *list , const size_t NMOM )
   @brief computes the momenta
   @warning only 2sin(q/2) momentum definition used at the moment
 */
void
compute_p_psq( double **p ,
	       double *psq ,
	       const struct veclist *list ,
	       const size_t NMOM ) ;

/**
   @fn void compute_WI( const struct PIdata *data , const double **p , const struct veclist *list , const size_t NMOM )
   @brief print to stdout the momentum space WIs
 */
void
compute_WI( const struct PIdata *data ,
	    const double **p ,
	    const struct veclist *list ,
	    const size_t NMOM) ;

/**
   @fn void correct_WI( struct PIdata *data , const correction_dir corr_dir , const struct veclist *list , const size_t NMOM )
   @brief perform WI correction on idx
   @warning overwrites data
 */
void
correct_WI( struct PIdata *data ,
	    const correction_dir corr_dir ,
	    const struct veclist *list ,
	    const size_t NMOM ) ;

/**
   @fn void WI_configspace_bwd( const struct PIdata *data , const struct site *lat )
   @brief configuration space WI test
   
   \f[ \frac{1}{V}\sum_{x,\mu,\nu} V_\mu(x) ||V_\nu(0) - V_\mu( x - \mu ) V_\nu( 0 )|| \f]
 */
void
WI_configspace_bwd( const struct PIdata *data ,
		    const struct site *lat ) ;


/**
   @fn void WI_configspace_fwd( const struct PIdata *data , const struct site *lat )
   @brief configuration space WI test
   
   \f[ \frac{1}{V}\sum_{x,\mu,\nu} V_\mu(x) ||V_\nu(0) - V_\mu( x - \mu ) V_\nu( 0 )|| \f]
 */
void
WI_configspace_fwd( const struct PIdata *data ,
		    const struct site *lat ) ;

/**
   @fn void WI_configspace_sym( const struct PIdata *data , const struct site *lat )
   @brief symmetric difference
 */
void
WI_configspace_sym( const struct PIdata *data ,
		    const struct site *lat ) ;

#endif
