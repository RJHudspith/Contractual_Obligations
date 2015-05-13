/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
*/
#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn const struct gamma CGmu( const struct gamma GAMMA_MU , const struct gamma *GAMMAS )
   @brief computes \f$ \gamma_t \gamma_y \f$
 */
struct gamma
CGmu( const struct gamma GAMMA_MU , 
      const struct gamma *GAMMAS ) ;

/**
   @fn const struct gamma CGmuT( const struct gamma Cgmu , const struct gamma *GAMMAS )
   @brief adjoint of CGmu \f$ = \gamma_t \left( C \gamma_\mu \right) \gamma_t \f$
 */
struct gamma
CGmuT( const struct gamma Cgmu , 
       const struct gamma *GAMMAS ) ;

/**
   @fn void gamma_mmul( struct gamma *__restrict a , const struct gamma b , const struct gamma c )
   @brief gamma multiplication
 */
void
gamma_mmul( struct gamma *__restrict a ,
	    const struct gamma b ,
	    const struct gamma c ) ;

/**
   @fn int make_gammas( struct gamma *GAMMA , const proptype prop )
   @brief defines the gamma matrices in whatever compiled convention
   @return #SUCCESS or #FAILURE
 */
int
make_gammas( struct gamma *__restrict GAMMA ,
	     const proptype prop ) ;

#endif





