/**
   @file contract_O2O2.h
   @brief prototype declarations for Baryon-Meson contractions
 */
#ifndef CONTRACT_O2O2_H
#define CONTRACT_O2O2_H

/**
   @fn void contract_O2O2( struct spinmatrix *P , const struct spinor L , const struct spinor S , const struct spinor B , const struct gamma *GAMMAS )
   @brief contract the Baryon-Meson operator
 */
void
contract_O2O2( struct spinmatrix *P ,
	       const struct spinor L ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS ) ;

#endif
