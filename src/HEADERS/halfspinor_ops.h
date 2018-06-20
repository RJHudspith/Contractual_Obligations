/**
   @file halfspinor_ops.h
   @brief prototype declarations for operations on halfspinor type
 */
#ifndef HALFSPINOR_OPS_H
#define HALFSPINOR_OPS_H

/**
   @fn void add_halfspinor( struct halfspinor *a , const struct halfspinor b )
   @brief atomically add two halfspinors s.t a = a + b
 */
void
add_halfspinor( struct halfspinor *a ,
		const struct halfspinor b ) ;

/**
   @fn void colormatrix_halfspinor( struct halfspinor *a , const double complex b[ NCNC ] , struct halfspinor c )
   @brief multiply a halfspinor on the left by a colormatrix
 */
void
colormatrix_halfspinor( struct halfspinor *a ,
			const double complex b[ NCNC ] ,
			struct halfspinor c ) ;

/**
   @fn void colormatrixdag_halfspinor( struct halfspinor *a , const double complex b[ NCNC ] , struct halfspinor c )
   @brief multiply a halfspinor on the left by the dagger of a color matrix
 */
void
colormatrixdag_halfspinor( struct halfspinor *a ,
			   const double complex b[ NCNC ] ,
			   struct halfspinor c ) ;

void
Fmunu_halfspinor( struct halfspinor *a ,
		  const double complex *b ,
		  const struct halfspinor c ) ;

/**
   @fn void halfspinor_iSaxpy( struct halfspinor *H , const struct halfspinor *S , const double fac )
   @brief computes the operation H += fac * ( I * S )
 */
void
halfspinor_iSaxpy( struct halfspinor *H ,
		   const struct halfspinor S ,
		   const double fac ) ;

/**
   @fn void halfspinor_multiply( struct halfspinor *a , const struct halfspinor b , const struct halfspinor c )
   @brief full spin-color multiply of two halfspinors a = b * c
 */
void
halfspinor_multiply( struct halfspinor *a ,
		     const struct halfspinor b ,
		     const struct halfspinor c ) ;

/**
   @fn void halfspinor_Saxpy( struct halfspinor *H , const struct halfspinor *S , const double fac )
   @brief computes the operation H += fac * S over
 */
void
halfspinor_Saxpy( struct halfspinor *H ,
		  const struct halfspinor S ,
		  const double fac ) ;

/**
   @fn void halfspinor_sigma_Saxpy( struct halfspinor *H , const struct halfspinor *S , const uint8_t imap[ NS ] , const size_t sigma_map[ NS ] ) 
   @brief Saxpy with multiplies by +/-1 or +/- I with implicit index shuffling in sigma map 
 */
void
halfspinor_sigma_Saxpy( struct halfspinor *H ,
			const struct halfspinor S ,
			const uint8_t sigma_map[ NS ] ,
			const uint8_t imap[ NS ] ) ;

/**
   @fn void zero_halfspinor( struct halfspinor *S )
   @brief zero a halfspinor
 */
void
zero_halfspinor( struct halfspinor *S ) ;

#endif
