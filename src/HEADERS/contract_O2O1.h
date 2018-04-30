/**
   @file contract_O2O1.h
   @brief prototype declaration for O2O1 pentaquark contraction
 */
#ifndef CONTRACT_O2O1_H
#define CONTRACT_O2O1_H

/**
   @fn void contract_O2O1( struct spinmatrix *P , const struct spinor U , const struct spinor D , const struct spinor S , const struct spinor B , const struct gamma OP1 , const struct gamma OP2 , const struct gamma *GAMMAS )
   @brief contraction of diquarks - Baryon-meson operator
 */
void
contract_O2O1( struct spinmatrix *P ,
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
