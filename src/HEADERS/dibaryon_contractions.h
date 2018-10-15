/**
   @file dibaryon_contractions.h
   @brief prototype declarations for SU(2) dibaryon contractions
 */
#ifndef DIBARYON_CONTRACTIONS_H
#define DIBARYON_CONTRACTIONS_H

/**
   @fn double complex dibaryon_contract( struct spinor S , const struct gamma *GAMMAS , const size_t GSRC , const size_t GSNK )
   @brief contract the SU(2) dibaryon
 */
double complex
dibaryon_contract( struct spinor S ,
		   const struct gamma *GAMMAS ,
		   const size_t GSRC ,
		   const size_t GSNK ) ;

// slower explicit version of the contraction
#if 0
double complex
dibaryon_contract2( struct spinor S ,
		   const struct gamma *GAMMAS ,
		    const size_t GSRC ,
		    const size_t GSNK ) ;
#endif

#endif
