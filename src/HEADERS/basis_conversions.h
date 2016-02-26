/**
   @file basis_conversions.h
   @brief prototype functions for converting between bases

   change from chiral into nrel basis
    
   chiral \f$ S = T * S * T^\dagger \f$
   here T:
   
           1/s    0       1/s     0
           0      1/s     0       1/s
           -1/s   0       1/s     0
           0      -1/s    0       1/s
    
   with \f$s = \sqrt{2} \f$
 */
#ifndef BASIS_CONVERSIONS_H
#define BASIS_CONVERSIONS_H

/**
   @fn void chiral_to_nrel( struct spinor *S )
   @brief convert a chiral spinor to a non-relativistic one by gamma basis rotation
 */
void
chiral_to_nrel( struct spinor *S ) ;

/**
   @fn void rotate_offdiag( struct spinor **S , const struct propagator *prop , const size_t Nprops )
   @brief rotate props if some are chiral and some are non-relativistic
 */
void 
rotate_offdiag( struct spinor **S ,
		const struct propagator *prop ,
		const size_t Nprops ) ;

#endif
