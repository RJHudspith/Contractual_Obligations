#ifndef PENTA_CONTRACTIONS_H
#define PENTA_CONTRACTIONS_H

/**
   @fn int pentas( double complex *result , const struct spinor L , const struct spinor S , const struct spinor bwdH , const struct gamma *GAMMAS , const size_t mu )
   @brief pentaquark contraction code for a udusb pentaquark
   @param L :: light quark propagator assumes u-d degeneracy
   @param S :: strange quark propagator
   @param bwdH :: backward-propagating heavy quark
   @param GAMMAS :: gamma matrices
   @param GSRC :: source gamma
   @param GSNK :: sink gamma
 */
int
pentas( double complex *result ,
	const struct spinor L , 
	const struct spinor S ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS ,
	const size_t GSRC ,
	const size_t GSNK ) ;

#endif
