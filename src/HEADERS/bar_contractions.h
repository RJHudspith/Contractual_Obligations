/**
   @file bar_contractions.h
   @brief prototype functions for baryon contractions
 */

#ifndef BAR_CONTRACTIONS_H
#define BAR_CONTRACTIONS_H

/**
   @fn const double complex baryon_contract( const struct spinor DiQ , const struct spinor S , const int d0 , const int d1 , const int d2 , const int d3 )
   @brief baryon contractions between Diquark and propagator S
 */
const double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const int d0 ,
		 const int d1 ,
		 const int d2 ,
		 const int d3 ) ;

/**
   @fn void cross_color_trace( struct spinor *__restrict DiQ , const struct spinor S )
   @brief cross-color trace of two propagators
   @warning overwrites DiQ with the cross-color trace expects DiQ to be \f$ Cg_\mu^{T} S Cg_\mu
 */
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S ) ;
#endif
