/**
   @file PImunu_projections.h
   @brief prototype functions for momentum-space projections
 */

#ifndef PIMUNU_PROJECTIONS_H
#define PIMUNU_PROJECTIONS_H

/**
   @fn void subtract_zeromom( struct PIdata *__restrict data , const int NMOM[1] )
   @brief subtracts the sum of the correlator from the data
 */
void
subtract_zeromom( struct PIdata *__restrict data ,
		  const int NMOM[1] ) ;

/**
   @fn void projection_fullmom( const struct PIdata *data , const double **p , const double *psq , const struct veclist *list , const int *__restrict NMOM , const char *outfile )
   @brief computes the transverse and longitudinal and their sum by projection

   \f[ \Pi^{1}(q^2) = \frac{1}{(ND-1)q^2}\left( \delta_{\mu\nu} - \frac{q_\mu q_\nu }{q^2}\right)\Pi_{\mu\nu}(q) \f]
   \f[ \Pi^{1}(q^2) = -\frac{q_\mu q_\nu}{q^4}\Pi_{\mu\nu}(q) \f]
 */
void
projection_fullmom( const struct PIdata *data ,
		    const double **p ,
		    const double *psq ,
		    const struct veclist *list ,
		    const int *__restrict NMOM ,
		    const char *outfile ) ;

/**
   @fn void projection_tonly( const struct PIdata *data , const double **p , const double *psq , const struct veclist *list , const int *__restrict NMOM , const char *outfile )
   @brief (0,0,0,pt) case is a special one
 */
void
projection_tonly( const struct PIdata *data ,
		  const double **p ,
		  const double *psq ,
		  const struct veclist *list ,
		  const int *__restrict NMOM ,
		  const char *outfile ) ;

#endif
