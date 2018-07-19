/**
   @file halfspinor_ops.h
   @brief prototype declarations for operations on halfspinor type
 */
#ifndef HALFSPINOR_OPS_SSE_H
#define HALFSPINOR_OPS_SSE_H

/**
   @fn void add_halfspinor( struct halfspinor *a , const struct halfspinor b )
   @brief atomically add two halfspinors s.t a = a + b
 */
void
add_halfspinor( struct halfspinor *a ,
		const struct halfspinor b ) ;

/**
   @fn void colormatrix_halfspinor( __m128d *pA , const __m128d *pB , const __m128d *pC )
   @brief multiply a halfspinor on the left by a colormatrix
 */
void
colormatrix_halfspinor( __m128d *pA ,
			const __m128d *pB ,
			const __m128d *pC ) ;

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
   @fn void sigmaB_halfspinor( struct halfspinor *S1 , const struct field Fmunu , const struct halfspinor S )
   @brief computes \sigma.B S
 */
void
sigmaB_halfspinor( struct halfspinor *S1 ,
		   const struct field Fmunu ,
		   const struct halfspinor S ) ;

/**
   @fn void zero_halfspinor( struct halfspinor *S )
   @brief zero a halfspinor
 */
void
zero_halfspinor( struct halfspinor *S ) ;

#endif
