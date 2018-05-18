/**
   @file grad_4.h
   @brief prototype declarations for higher order derivatives
 */
#ifndef GRAD_4_H
#define GRAD_4_H

/**
   @fn void grad4( struct halfspinor *der , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const size_t mu )
   @brief computes the fourth-order derivative of S at site i and timeslice t
 */
void
grad4( struct halfspinor *der ,
       const struct halfspinor *S ,
       const struct field *Fmunu ,
       const size_t i ,
       const size_t t ,
       const size_t mu ) ;

/**
   @fn void grad_sqsq( struct halfspinor *der2 , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t )
   @brief computes grad^2_nu grad^2_mu summed over mu and nu at site i and timeslice t
 */
void
grad_sqsq( struct halfspinor *der2 ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t ) ;

#endif
