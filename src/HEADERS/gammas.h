/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
 */

#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn int gamma_matrix( double complex *g, int *ig , const int mu )
   @brief reshuffles indices (implicit gamma multiplication)
 */
int 
gamma_matrix( double complex *g, 
	      int *ig , 
	      const int mu ) ;

#endif





