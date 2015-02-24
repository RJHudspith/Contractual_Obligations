/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
*/
#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn void gamma_matrix( struct gamma *G , const proptype prop )
   @brief defines the gamma matrices in whatever compiled convention
 */
void
make_gammas( struct gamma *GAMMA ,
	     const proptype prop ) ;

#endif





