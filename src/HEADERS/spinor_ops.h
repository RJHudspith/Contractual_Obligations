/**
   @file spinor_ops.h
   @brief prototype functions for various spinor operations
 */
#ifndef SPINOR_OPS_H
#define SPINOR_OPS_H

#include "../config.h"

#ifdef HAVE_EMMINTRIN_H

// we include the SSE instructions where possible
#include "spinor_ops_SSE.h"

#else

/**
   @fn void add_spinors( struct spinor *A , const struct spinor B )
   @brief atomically add two spinors A += B
 */
void
add_spinors( struct spinor *A ,
	     const struct spinor B ) ;


/**
   @fn void colortrace_spinor( void *S1 , const void *S2 )
   @brief trace the color indices of our spinor
 */
void
colortrace_spinor( void *S1 ,
		   const void *S2 ) ;

/**
   @fn void equate_spinor( void *S1 , const void *S2 )
   @brief equate two spinors
 */
void
equate_spinor( void *S1 ,
	       const void *S2 ) ;

/**
   @fn void equate_spinor_minus( void *S1 , const void *S2 )
   @brief equate one spinor to the minus of another
 */
void
equate_spinor_minus( void *S1 ,
		     const void *S2 ) ;

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
   @fn void identity_spinor( struct spinor *__restrict res )
   @brief set a spinor to be deltas in both color and spin
 */
void
identity_spinor( struct spinor *__restrict res ) ;

/**
   @fn void spinor_gauge( struct spinor *__restrict res , const struct spinor S , const double complex link[ NCNC ] ) 
   @brief multiplies a spinor with a link matrix \f$ res = S * link \f$
 */
void
spinor_gauge( struct spinor *__restrict res ,
	      const struct spinor S ,
	      const double complex link[ NCNC ] ) ;

/**
   @fn void spinmul_atomic_left( struct spinor *A , const struct spinor B )
   @brief atomically left multiply spinor A by spinor B ( A = B * A )
 */
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B ) ;

/**
   @fn void spinmul_atomic_left( struct spinor *A , const struct spinor B )
   @brief atomically left multiply spinor A by spinor B ( A = A * B )
 */
void
spinmul_atomic_right( struct spinor *A ,
		      const struct spinor B ) ;

/**
   @fn void spinor_zero( void *S )
   @brief zero a spinor over the timeslice
   @warning is threaded
 */
void
spinor_zero( void *S ) ;

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
   @fn void spinor_Saxpy( struct spinor *A , const double S , const struct spinor *B )
   @brief compute A = A + S*B
*/
void
spinor_Saxpy( struct spinor *A ,
	      const double S ,
	      const struct spinor *B ) ;

/**
   @fn void spinor_zero_site( void *S )
   @brief zero a spinor on a site
 */
void
spinor_zero_site( void *S ) ;

/**
   @fn void spintrace( void *S , const void *S2 )
   @brief traces over spin indices of S2 into an #NCNC flattened color matrix S
   @param S :: color matrix
   @param S2 :: spinor
 */
void
spintrace( void *S ,
	   const void *S2 ) ;

/**
   @fn void sub_spinors( struct spinor *A , const struct spinor B )
   @brief atomically subtract two spinors A -= B
 */
void
sub_spinors( struct spinor *A ,
	     const struct spinor B ) ;

/**
   @fn void sumprop( void* SUM , const void* S )
   @brief sum a propagator over a timeslice
 */
void
sumprop( void *SUM ,
	 const void *S ) ;

/**
   @fn void sumwalls( struct spinor *SUM , const struct spinor **S , const size_t Nprops )
   @brief sum over Nprops of spinors over #LCU
   @warning is threaded
 */
void
sumwalls( struct spinor *SUM ,
	  const struct spinor **S ,
	  const size_t Nprops ) ;

/**
   @fn struct spinor transpose_spinor( const struct spinor S )
   @brief transpose a spinor @S
   #return the transpose of S
 */
struct spinor
transpose_spinor( const struct spinor S ) ;

#endif

#endif
