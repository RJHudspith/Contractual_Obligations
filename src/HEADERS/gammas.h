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
   @fn void gamma_mul_minus1( struct gamma *G )
   @brief multiply a gamma matrix by -1
 */
void
gamma_mul_minus1( struct gamma *G ) ;

/**
   @fn void gamma_mul_minusi( struct gamma *G )
   @brief multiply a gamma matrix by -i
 */
void
gamma_mul_minusi( struct gamma *G ) ;

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
   @fn struct gamma gt_Gconj_gt( const struct gamma G , const struct gamma gt )
   @brief computes \f$ \gamma_t G^* \gamma_t \f$
 */
struct gamma
gt_Gconj_gt( const struct gamma G , 
	     const struct gamma gt ) ;

/**
   @fn struct gamma gt_Gdag_gt( const struct gamma G , const struct gamma gt )
   @brief computes \f$ \gamma_t G^\dagger \gamma_t \f$
 */
struct gamma
gt_Gdag_gt( const struct gamma G , 
	    const struct gamma gt ) ;

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

/**
   @fn int setup_gamma_2( struct gamma *GAMMA , const proptype basis1 , const proptype basis2 )
   @brief setup the gammas dependent on the basis of prop1 or prop2
   @return #SUCCESS or #FAILURE
 */
int
setup_gamma_2( struct gamma *GAMMAS ,
	       const proptype basis1 ,
	       const proptype basis2 ) ;

/**
   @fn int setup_gamma_3( struct gamma *GAMMAS , const proptype basis1 , const proptype basis2 , const proptype basis3 )
   @brief setup the gammas dependent on the basis of prop1 or prop2 or prop3
   @return #SUCCESS or #FAILURE
 */
int
setup_gamma_3( struct gamma *GAMMAS ,
	       const proptype basis1 ,
	       const proptype basis2 , 
	       const proptype basis3 ) ;

#endif





