/**
   @file diquark_contraction.h
   @brief prototype declarations for diquark contraction
 */
#ifndef DIQUARK_CONTRACTION_H
#define DIQUARK_CONTRACTION_H

/**
   @fn double complex diquark( struct spinor S1 , struct spinor S2 , const struct gamma C_GSRC , const struct gamma C_GSNK ) 
   @brief contract two diquarks
 */
double complex
diquark( struct spinor S1 , 
	 struct spinor S2 ,
	 const struct gamma C_GSRC , 
	 const struct gamma C_GSNK ,
	 const struct gamma G5 ) ;

#endif
