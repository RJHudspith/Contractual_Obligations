/**
   @file basis_conversions.h
   @brief prototype functions for converting between bases
 */
#ifndef BASIS_CONVERSIONS_H
#define BASIS_CONVERSIONS_H

/**
   @fn void chiral_to_nrel( struct spinor *S )
   @brief convert chiral basis to non-relativistic

   change from chiral into nrel basis
    
   chiral \f$ S = T * S * T^\dagger \f$
   here T:
   
           1/s    0       1/s     0
           0      1/s     0       1/s
           -1/s   0       1/s     0
           0      -1/s    0       1/s
    
   with \f$s = \sqrt{2} \f$
 */
void
chiral_to_nrel( struct spinor *S ) ;

/**
   @fn void nrel_rotate_slice( struct spinor *S )
   @brief rotate a time-slice worth of spinor S to NREL basis
 */
void
nrel_rotate_slice( struct spinor *S ) ;

/**
   @fn void rotate_offdiag_2( struct spinor *S1 , const proptype basis1 , struct spinor *S2 , const proptype basis2 )
   @brief rotate either S1 or S2 to the non-rel basis if one is chiral and the other non-relativistic
 */
void
rotate_offdiag_2( struct spinor *S1 ,
		  const proptype basis1 ,
		  struct spinor *S2 ,
		  const proptype basis2 ) ;

/**
   @fn void rotate_offdiag_3( struct spinor *S1 , const proptype basis1 , struct spinor *S2 , const proptype basis2 , struct spinor *S3 , const proptype basis3 )
   @brief rotate chiral to nrel basis for any of the three props given
 */
void
rotate_offdiag_3( struct spinor *S1 ,
		  const proptype basis1 ,
		  struct spinor *S2 ,
		  const proptype basis2 , 
		  struct spinor *S3 ,
		  const proptype basis3 ) ;

#endif
