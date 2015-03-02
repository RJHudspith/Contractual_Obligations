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
   @fn void gauge_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
gauge_spinor( struct spinor *res ,  
	      const double complex link[ NCNC ] ,
	      const struct spinor S ) ;

/**
   @fn void gauge_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
spinor_gauge( struct spinor *res ,
	      const struct spinor S ,
	      const double complex link[ NCNC ] ) ;

/**
   @fn void gaugedag_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a daggered link res = link^{\dagger} S
 */
void
gaugedag_spinor( struct spinor *res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

void
spinordag_gauge( struct spinor *res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

void
gauge_spinordag( struct spinor *res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

void
spinor_gaugedag( struct spinor *res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

/**
   @brief trace of the product of two matrices Tr[ A*B ]
 */
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B ) ;

/**
   @brief atomic left multiply a spinor by a gamma matrix
 */
void
gamma_mul_l( struct spinor *res ,
	     const struct gamma GAMMA ) ;

/**
   @brief atomic right multiply a spinor by a gamma matrix
 */
void
gamma_mul_r( struct spinor *res ,
	     const struct gamma GAMMA ) ;

/**
   @brief conjugate transpose a spinor
 */
void
adjoint_spinor( struct spinor *adj ,
		const struct spinor S ) ;

/**
   @brief computes gamma_5 adj( S ) gamma_5, puts result in adj
 */
void
full_adj( struct spinor *adj ,
	  const struct spinor S ,
	  const struct gamma G5 ) ;

#endif
