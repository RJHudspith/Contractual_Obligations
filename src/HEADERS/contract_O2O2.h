/**
   @file contract_O2O2.h
   @brief prototype declarations for Baryon-Meson contractions
 */
#ifndef CONTRACT_O2O2_H
#define CONTRACT_O2O2_H

/**
   @fn void contract_O2O2( struct spinmatrix *P , const struct spinor U , const struct spinor D , const struct spinor S , const struct spinor B , const struct gamma OP1 , const struct gamma OP2 , const struct gamma *GAMMAS )
   @brief contract the Baryon-Meson operator
 */
void
contract_O2O2( struct spinmatrix *P ,
	       double complex **F ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS ,
	       const uint8_t **loc ) ;

#endif
