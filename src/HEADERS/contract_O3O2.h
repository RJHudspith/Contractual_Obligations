/**
   @file contract_O2O3.h
   @brief prototype declarations for the first and second Baryon-meson ops
 */
#ifndef CONTRACT_O3O2_H
#define CONTRACT_O3O2_H

/**
   @fn void contract_O3O2( struct spinmatrix *P , const struct spinor U , const struct spinor D , const struct spinor S , const struct spinor B , const struct gamma OP1 , const struct gamma OP2 , const struct gamma *GAMMAS )
   @brief baryon-meson mixing contraction
 */
void
contract_O3O2( struct spinmatrix *P ,
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
