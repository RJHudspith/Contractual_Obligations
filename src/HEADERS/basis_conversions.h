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

#endif
