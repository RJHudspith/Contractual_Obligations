/**
   @file correlators.h
   @brief correlation function calculation prototypes
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H


/**
   @fn double complex local_meson_correlator( const struct spinor S1 , const struct spinor S2 , const int *sadj , const int *iadj, const int mu1, const int mu2 )
   @return C(t) from spin-color trace
 */
double complex 
local_meson_correlator( const struct spinor S1 , 
			const struct spinor S2 , 
			const int *sadj ,
			const int *iadj , 
			const int mu1 , 
			const int mu2 ) ;

/**
   @fn double complex local_conserved_meson_correlator( const struct spinor S1 , const struct spinor S2 , const int *sadj , const int *iadj, const int mu1, const int mu2 )
   @brief local-conserved correlator
 */
double complex 
local_conserved_meson_correlator( const struct spinor S1 , 
				  const struct spinor S2 , 
				  const int *sadj ,
				  const int *iadj , 
				  const int mu1 , 
				  const int mu2 ) ;


#endif
