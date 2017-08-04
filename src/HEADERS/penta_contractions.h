/**
   @file penta_contractions.h
   @brief prototype declarations for pentaquark contractions
 */
#ifndef PENTA_CONTRACTIONS_H
#define PENTA_CONTRACTIONS_H

/**
   @fn size_t idx( const size_t b , const size_t bp , const size_t c , const size_t cp , const size_t g , const size_t gp , const size_t h , const size_t hp )
   @brief return a linearised index from a tensor with 8 indices
 */
size_t
idx( const size_t b , const size_t bp ,
     const size_t c , const size_t cp ,
     const size_t g , const size_t gp ,
     const size_t h , const size_t hp ) ;

/**
   @fn size_t idx( const size_t b , const size_t bp , const size_t c , const size_t cp , const size_t g , const size_t gp , const size_t h , const size_t hp )
   @brief return a linearised index from a tensor with 6 indices
 */
size_t
idx2( const size_t b , const size_t bp ,
      const size_t c , const size_t cp ,
      const size_t g , const size_t gp ) ;

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
