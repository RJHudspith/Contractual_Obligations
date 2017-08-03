/**
   @file penta_contractions.h
   @brief prototype declarations for pentaquark contractions
 */
#ifndef PENTA_CONTRACTIONS_H
#define PENTA_CONTRACTIONS_H

/**
   @fn size_t idx( const size_t a , const size_t ap , const size_t b , const size_t bp , const size_t c , const size_t cp , const size_t d , const size_t dp )
   @brief return a linearised index from a tensor with 8 indices
 */
size_t
idx( const size_t a , const size_t ap ,
     const size_t b , const size_t bp ,
     const size_t c , const size_t cp ,
     const size_t d , const size_t dp ) ;

/**
   @fn int pentas( double complex *result , const struct spinor L , const struct spinor S , const struct spinor bwdH , const struct gamma *GAMMAS )
   @brief pentaquark contraction code for a udusb pentaquark
   @param L :: light quark propagator assumes u-d degeneracy
   @param S :: strange quark propagator
   @param bwdH :: backward-propagating heavy quark
   @param GAMMAS :: gamma matrices
 */
int
pentas( double complex *result ,
	const struct spinor L , 
	const struct spinor S ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS ) ;

#endif
