/**
   @file correlators.h
   @brief correlation function calculation prototypes
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn double complex pion_correlator( const struct spinor S1 , const struct spinor S2 ) 
   @brief computes the Pion correlation function for a site
   @param S1 :: propagator 1
   @param S2 :: propagator 2
   @return the spin-color trace
 */
double complex 
pion_correlator( const struct spinor S1 ,
		 const struct spinor S2 ) ;

/**
   @fn double complex local_meson_correlator( const struct spinor S1 , const struct spinor S2 , const double complex *sadj , const int *iadj, const int mu1, const int mu2 )
   @return C(t) from spin-color trace
 */
double complex 
local_meson_correlator( const struct spinor S1 , 
			const struct spinor S2 , 
			const double complex *sadj ,
			const int *iadj , 
			const int mu1 , 
			const int mu2 ) ;

/**
   @fn double complex local_conserved_meson_correlator( const struct spinor S1 , const struct spinor S2 , const double complex *sadj , const int *iadj, const int mu1, const int mu2 )
   @brief local-conserved correlator
 */
double complex 
local_conserved_meson_correlator( const struct spinor S1 , 
				  const struct spinor S2 , 
				  const double complex *sadj ,
				  const int *iadj , 
				  const int mu1 , 
				  const int mu2 ) ;
#endif
