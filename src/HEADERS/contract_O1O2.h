/**
   @file contract_O1O2.h
   @brief prototype declarations for O1O2 pentaquark contractions
 */
#ifndef CONTRACT_O1O2_H
#define CONTRACT_O1O2_H

/**
   @fn void contract_O1O2( struct spinmatrix *P , const struct spinor L , const struct spinor S , const struct spinor B , const struct gamma *GAMMAS )
   @brief contract the diquarks - Baryon-Meson operators
 */
void
contract_O1O2( struct spinmatrix *P ,
	       const struct spinor L ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS ) ;

#endif
