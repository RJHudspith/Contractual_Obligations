/**
   @file spinor_ops_SSE.h
   @brief intrinsics spinor operations
 */

#ifndef SPINOR_OPS_SSE_H
#define SPINOR_OPS_SSE_H

#ifdef HAVE_EMMINTRIN_H

/**
   @fn void equate_spinor( void *S1 , const void *S2 )
   @brief equate two spinors
 */
void
equate_spinor( void *S ,
	       const void *S2 ) ;

/**
   @fn void equate_spinor_minus( void *S1 , const void *S2 )
   @brief equate one spinor to the minus of another
 */
void
equate_spinor_minus( void *mS ,
		     const void *S ) ;

/**
   @fn void flipsign_spinor( void *S ) 
   @brief flips the sign of the spinor s
 */
void
flipsign_spinor( void *S ) ;

/**
   @fn void gauge_spinor( struct spinor *__restrict res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
gauge_spinor( struct spinor *__restrict res ,
	      const double complex link[ NCNC ] ,
	      const struct spinor S ) ;

/**
   @fn void gaugedag_spinor( struct spinor *__restrict res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a daggered link \f$ res = link^{\dagger} S \f$
 */
void
gaugedag_spinor( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

/**
   @fn void gauge_spinordag( struct spinor *__restrict res , const double complex link[ NCNC ] , const struct spinor S ) 
   @brief computes \f$ res = U.S^{\dagger} \f$
 */
void
gauge_spinordag( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

/**
   @fn void spinor_gauge( struct spinor *__restrict res , const struct spinor S , const double complex link[ NCNC ] ) 
   @brief multiplies a spinor with a link matrix \f$ res = S * link \f$
 */
void
spinor_gauge( struct spinor *__restrict res ,  
	      const struct spinor S ,
	      const double complex link[ NCNC ] ) ;

/**
   @fn void spinordag_gauge( struct spinor *__restrict res , const struct spinor S , const double complex link[ NCNC ] )
   @brief computes \f$ res = S^{\dagger}.U \f$
 */
void
spinordag_gauge( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

/**
   @fn void spinor_gaugedag( struct spinor *__restrict res , const struct spinor S , const double complex link[ NCNC ] )
   @brief computes \f$ res = S.U^{\dagger} \f$
 */
void
spinor_gaugedag( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

/**
   @fn void sumprop( void* SUM , const void* S )
   @brief sum a propagator over a timeslice
 */
void
sumprop( void *SUM ,
	 const void *S ) ;

/**
   @fn void spinmul_atomic_left( struct spinor *A , const struct spinor B )
   @brief atomically left multiply spinor A by spinor B ( A = B * A )
 */
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B ) ;

/**
   @fn void spinor_zero( void *S )
   @brief zero a spinor over the timeslice
   @warning is threaded
 */
void
spinor_zero( void *S ) ;

#endif

#endif
