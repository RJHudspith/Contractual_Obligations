/**
   @file tetra_contractions.h
   @brief prototype declarations for tetraquark contractions
 */
#ifndef TETRA_CONTRACTIONS_H
#define TETRA_CONTRACTIONS_H

/**
   @fn int tetras( double complex result[ TETRA_NOPS ] , const struct spinor L1 , const struct spinor L2 , const struct spinor bwdH , const struct gamma *GAMMAS , const size_t mu , const GLU_bool L1L2_degenerate )
   @brief perform all tetraquark contractions
   @return #SUCCES or #FAILURE
 */
int
tetras( double complex result[ TETRA_NOPS ] ,
	const struct spinor L1 , 
	const struct spinor L2 ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS , 
	const size_t mu , 
	const GLU_bool L1L2_degenerate ) ;

#endif
