/**
   @file matrix_ops.h
   @brief prototype functions for general matrix operations
 */
#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication a = b * c
 */
void 
multab( double complex a[ NCNC ] , 
	const double complex b[ NCNC ] , 
	const double complex c[ NCNC ] ) ;

/**
   @fn void multabdag( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication a = b^{\dagger} * c
 */
void 
multabdag( double complex a[ NCNC ] , 
	   const double complex b[ NCNC ] , 
	   const double complex c[ NCNC ] ) ;

/**
   @fn void gauge_spinor( struct spinor *res , const struct spinor S , const double complex link[ NCNC ] )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
gauge_spinor( struct spinor *res ,  
	      const struct spinor S ,
	      const double complex link[ NCNC ] ) ;

/**
   @fn void gaugedag_spinor( struct spinor *res , const struct spinor S , const double complex link[ NCNC ] )
   @brief multiplies a spinor with a daggered link res = link^{\dagger} S
 */
void
gaugedag_spinor( struct spinor *res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

#endif
