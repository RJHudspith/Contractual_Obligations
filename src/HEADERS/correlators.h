/**
   @file correlators.h
   @brief correlation function calculation prototypes
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn double complex local_meson_correlator_singlet( const struct spinor S1 , const struct gamma ADJ , const struct gamma SRC , const struct gamma SNK )
   @param S1 :: propagator solution at a site

   @return spin-color trace
 */
double complex 
local_meson_correlator_singlet( const struct spinor S1 , 
				const struct gamma ADJ ,
				const struct gamma SRC ,
				const struct gamma SNK ) ;

/**
   @fn double complex local_meson_correlator( const struct spinor S1 , const struct spinor S2 , const struct gamma ADJ , const struct gamma SRC , const struct gamma SNK )
   @param S1 :: propagator solution at a site
   @param S2 :: second propagator solution at a site

   @return spin-color trace
 */
double complex 
local_meson_correlator( const struct spinor S1 , 
			const struct spinor S2 , 
			const struct gamma ADJ ,
			const struct gamma SRC ,
			const struct gamma SNK ) ;

#endif
