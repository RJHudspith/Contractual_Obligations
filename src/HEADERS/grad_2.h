/**
   @file grad_2.h
   @brief second order derivative prototype declarations
 */
#ifndef GRAD_2_H
#define GRAD_2_H

/**
   @fn void gradsq( struct halfspinor *der2 , const struct halfspinor *S , const size_t i , const size_t t )
   @brief computes grad^2 of S at site i on timeslice t
 */
void
gradsq( struct halfspinor *der2 ,
	const struct halfspinor *S ,
	const size_t i ,
	const size_t t ) ;

/**
   @fn void grad_sq_LCU( struct halfspinor *der2 , const struct halfspinor *S , const size_t t )
   @brief computes grad^2 over #LCU and puts into der2
 */
void
grad_sq_LCU( struct halfspinor *der2 ,
	     const struct halfspinor *S ,
	     const size_t t ) ;

/**
   @fn void gradsq_imp( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t )
   @brief computes the improved grad^2 at site i on timeslice t
 */
void
gradsq_imp( struct halfspinor *der ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ) ;

/**
   @fn void gradsq_imp_sigmaB( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t )
   @brief computes \grad^2 \sigma.B S
*/
void
gradsq_imp_sigmaB( struct halfspinor *der ,
		   const struct halfspinor *S ,
		   const struct field *Fmunu ,
		   const size_t i ,
		   const size_t t ) ;

/**
   @fn void sigmaB_gradsq_imp( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t )
   @brief compute \sigma.b \grad^2 S
*/
void
sigmaB_gradsq_imp( struct halfspinor *der ,
		   const struct halfspinor *S ,
		   const struct field *Fmunu ,
		   const size_t i ,
		   const size_t t ) ;

#endif
