/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
*/
#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn void gamma_mmul( struct gamma *a , const struct gamma b , const struct gamma c )
   @brief gamma multiplication
 */
void
gamma_mmul( struct gamma *a ,
	    const struct gamma b ,
	    const struct gamma c ) ;

/**
   @fn int gamma_matrix( struct gamma *G , const proptype prop )
   @brief defines the gamma matrices in whatever compiled convention
   @return #SUCCESS or #FAILURE
 */
int
make_gammas( struct gamma *GAMMA ,
	     const proptype prop ) ;

#endif





