/**
   @file tetra_contractions.h
   @brief prototype declarations for tetraquark contractions
 */
#ifndef TETRA_CONTRACTIONS_H
#define TETRA_CONTRACTIONS_H

/**
   @fn size_t element( const size_t a , const size_t b , const size_t c , const size_t d ) 
   @brief given certain color components returns the linearised index
 */
size_t
element( const size_t a , 
	 const size_t b , 
	 const size_t c , 
	 const size_t d ) ;

/**
   @fn void get_abcd( size_t *a , size_t *b , size_t *c , size_t *d , const size_t abcd )
   @brief given the linearised index @abcd give back the correct color components 
 */
void
get_abcd( size_t *a , 
	  size_t *b , 
	  size_t *c , 
	  size_t *d , 
	  const size_t abcd ) ;

/**
   @fn void precompute_block( struct block *C1 , const struct spinor S1 , const struct gamma G1 , const struct spinor S2 , const struct gamma G2 )
   @brief precomputes spinmatrices with exposed abcd color indices
 */
void
precompute_block( struct block *C1 ,
		  const struct spinor S1 ,
		  const struct gamma G1 ,
		  const struct spinor S2 ,
		  const struct gamma G2 ) ;

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
