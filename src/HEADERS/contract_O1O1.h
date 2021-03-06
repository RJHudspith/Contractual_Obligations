/**
   @file contract_O1O1.h
   @brief diquarks contraction for the pentaquark
 */
#ifndef CONTRACT_O1O1_H
#define CONTRACT_O1O1_H
  
/**
   @fn void contract_O1O1( struct spinmatrix *P , const struct spinor U , const struct spinor D , const struct spinor S , const struct spinor bwdH , const struct gamma OP1 , const struct gamma OP2 , const struct gamma *GAMMAS )
   @brief contract the diquarks for the pentaquark
 */
void
contract_O1O1( struct spinmatrix *P ,
	       double complex **F ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor bwdH ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 , 
	       const struct gamma *GAMMAS ,
	       const uint8_t **loc ) ;

#endif
