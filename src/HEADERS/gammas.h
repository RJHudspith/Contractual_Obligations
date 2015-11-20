/**
   @file gammas.h
   @brief prototype functions for gamma matrix ops
*/
#ifndef GAMMAS_H
#define GAMMAS_H

/**
   @fn const struct gamma CGmu( const struct gamma G , const struct gamma *GAMMAS )
   @brief adjoint of CGmu \f$ = \gamma_t \left( C \gamma_\mu \right) \gamma_t \f$
   Where C is \f$C=i\gamma_y\gamma_t \f$ as in Gattringer & Lang Eq.A.24
 */
struct gamma
CGmu( const struct gamma G , 
       const struct gamma *GAMMAS ) ;

/**
   @fn struct gamma gamma_conj( const struct gamma G )
   @brief conjugate a gamma matrix
 */
struct gamma
gamma_conj( const struct gamma G ) ;

/**
   @fn struct gamma gamma_dag( const struct gamma G )
   @brief dagger a gamma matrix
 */
struct gamma
gamma_dag( const struct gamma G ) ;

/**
   @fn void gamma_mmul( struct gamma *__restrict a , const struct gamma b , const struct gamma c )
   @brief gamma multiplication
 */
void
gamma_mmul( struct gamma *__restrict a ,
	    const struct gamma b ,
	    const struct gamma c ) ;

/**
   @fn void gamma_muli( struct gamma *G )
   @brief multiply a gamma matrix by i
 */
void
gamma_muli( struct gamma *G ) ;

/**
   @fn struct gamma gamma_transpose( const struct gamma G )
   @brief transpose a gamma matrix
 */
struct gamma
gamma_transpose( const struct gamma G ) ;

/**
   @fn uint8_t gconj( const uint8_t )
   @brief return a conjugate element of our gammas
 */
uint8_t
gconj( const uint8_t ) ;

/**
   @fn struct gamma CGmu( const struct gamma G , const struct gamma *GAMMAS )
   @brief computes \f$ \gamma_t \gamma_y G \f$
 */
struct gamma
gt_Gconj_gt( const struct gamma G , 
	     const struct gamma *GAMMAS ) ;

/**
   @fn struct gamma gt_Gdag_gt( const struct gamma Cgmu , const struct gamma *GAMMAS )
   @brief computes \f$ \gamma_t G^\dagger \gamma_t \f$
 */
struct gamma
gt_Gdag_gt( const struct gamma G , 
	    const struct gamma *GAMMAS ) ;

/**
   @fn int make_gammas( struct gamma *GAMMA , const proptype prop )
   @brief defines the gamma matrices in whatever compiled convention
   @return #SUCCESS
 */
int
make_gammas( struct gamma *__restrict GAMMA ,
	     const proptype prop ) ;

/**
   @fn picture_gamma( const struct gamma G ) 
   @brief print to stdout gamma G
 */
void
picture_gamma( const struct gamma G ) ;

#endif





