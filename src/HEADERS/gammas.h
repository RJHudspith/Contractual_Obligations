/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
*/
#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn void gamma_matrix( struct gamma *G )
   @brief defines the gamma matrices in whatever compiled convention
 */
void
make_gammas( struct gamma *GAMMA ) ;

/**
   @fn void make_gammas_nrel( struct gamma *GAMMA )
   @brief gamma matrices in non-relativistic basis
 */
void
make_gammas_nrel( struct gamma *GAMMA ) ;

#endif





