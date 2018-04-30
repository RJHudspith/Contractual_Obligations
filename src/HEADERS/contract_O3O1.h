/**
   @file contract_O3O1.h
   @brief prototype declarations for the first and second Baryon-meson ops
 */
#ifndef CONTRACT_O3O1_H
#define CONTRACT_O3O1_H

/**
   @fn void contract_O3O1( struct spinmatrix *P , const struct spinor U , const struct spinor D , const struct spinor S , const struct spinor B , const struct gamma OP1 , const struct gamma OP2 , const struct gamma *GAMMAS )
   @brief baryon-meson mixing contraction
 */
void
contract_O3O1( struct spinmatrix *P ,
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
