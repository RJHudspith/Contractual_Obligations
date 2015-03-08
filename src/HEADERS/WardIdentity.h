/**
   @file WardIdentity.h
   @brief prototype functions for ward identity checking
 */

#ifndef WARDIDENTITY_H
#define WARDIDENTITY_H

/**
   @fn void compute_p_psq( double **p , double *psq , const struct veclist *list , const int NMOM )
   @brief computes the momenta
   @warning only 2sin(q/2) momentum definition used at the moment
 */
void
compute_p_psq( double **p ,
	       double *psq ,
	       const struct veclist *list ,
	       const int NMOM ) ;

/**
   @fn void compute_WI( const struct PIdata *data , const double **p , const struct veclist *list , const int NMOM )
   @brief print to stdout the momentum space WIs
 */
void
compute_WI( const struct PIdata *data ,
	    const double **p ,
	    const struct veclist *list ,
	    const int NMOM) ;

/**
   @fn void correct_WI( struct PIdata *data , const correction_dir corr_dir , const struct veclist *list , const int NMOM )
   @brief perform WI correction on idx
   @warning overwrites data
 */
void
correct_WI( struct PIdata *data ,
	    const correction_dir corr_dir ,
	    const struct veclist *list ,
	    const int NMOM ) ;

#endif
