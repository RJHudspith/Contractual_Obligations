/**
   @file Ospinor.h
   @brief opposite-ordered spinor operations prototype functions
 */
#ifndef OSPINOR_H
#define OSPINOR_H

/**
   @fn struct Ospinor spinor_to_Ospinor( const struct spinor S )
   @brief convert a spinor to an Ospinor
 */
struct Ospinor
spinor_to_Ospinor( const struct spinor S ) ;

/**
   @fn void gamma_mul_r_Ospinor( struct Ospinor *S , const struct gamma Gamma )
   @brief right multiply an Ospinor by a gamma matrix
 */
void
gamma_mul_r_Ospinor( struct Ospinor *S ,
		     const struct gamma Gamma ) ;

/**
   @fn void gamma_mul_r_OspinorT( struct Ospinor *S , const struct gamma Gamma )
   @brief right multiply an Ospinor by a gamma matrix and transposes the spin indices
 */
void
gamma_mul_r_OspinorT( struct Ospinor *S ,
		      const struct gamma Gamma ) ;

#endif
