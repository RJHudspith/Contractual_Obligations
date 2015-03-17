/**
   @file spinor_ops.h
   @brief prototype functions for various spinor operations
 */
#ifndef SPINOR_OPS_H
#define SPINOR_OPS_H

/**
   @fn void add_spinors( double complex *__restrict SUM , const double complex *__restrict S )
   @brief atomically add two spinors SUM += S
 */
void
add_spinors( double complex *__restrict SUM ,
	     const double complex *__restrict S ) ;

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
   @fn void sumprop( struct spinor *__restrict SUM , struct spinor *__restrict S )
   @brief sum a propagator over a timeslice
 */
void
sumprop( struct spinor *__restrict SUM ,
	 const struct spinor *__restrict S ) ;

/**
   @fn void spinmul_atomic_left( struct spinor *A , const struct spinor B )
   @brief atomically left multiply spinor A by spinor B ( A = B * A )
 */
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B ) ;

/**
   @fn void zero_spinor( double complex *__restrict S )
   @brief set all elements of S to 0
 */
void
zero_spinor( double complex *__restrict S ) ;

#endif
