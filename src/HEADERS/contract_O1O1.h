/**
   @file contract_O1O1.h
   @brief diquarks contraction for the pentaquark
 */
#ifndef CONTRACT_O1O1_H
#define CONTRACT_O1O1_H
  
/**
   @fn void contract_O1O1( struct spinmatrix *P , const struct spinor L , const struct spinor S , const struct spinor bwdH , const struct gamma *GAMMAS )
   @brief contract the diquarks for the pentaquark
 */
void
contract_O1O1( struct spinmatrix *P ,
	       const struct spinor L ,
	       const struct spinor S ,
	       const struct spinor bwdH ,
	       const struct gamma *GAMMAS ) ;

#endif
