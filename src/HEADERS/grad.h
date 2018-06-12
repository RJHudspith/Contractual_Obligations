/**
   @file grad.h
   @brief improved gradient computations
 */
#ifndef DERIVS_H
#define DERIVS_H

/**
   @fn void grad_imp_LCU( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const size_t t , const size_t mu )
   @brief compute the improved gradient of S over the whole timeslice and put into der
   @warning this code is never used but it was too cute to delete
 */
void
grad_imp_LCU( struct halfspinor *der ,
	      const struct halfspinor *S ,
	      const struct field *Fmunu ,
	      const size_t t ,
	      const size_t mu ) ;

/**
   @fn void FMUNU_grad_imp( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const double U_0 , const size_t i , const size_t t , const size_t mu , const size_t Fmunu_idx )
   @brief computes the improved derivative of S and left multiplies by Fmunu[i].O[Fmunu_idx] at site index i on timeslice t
 */
void
FMUNU_grad_imp( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx ) ;
/**
   @fn void grad_imp_FMUNU( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const double U_0 ,const size_t i , const size_t t , const size_t mu , const size_t Fmunu_idx ) 
   @brief computes the improved derivative of Fmunu.S for site i
 */
void
grad_imp_FMUNU( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx ) ;

#endif
