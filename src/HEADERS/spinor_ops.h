/**
   @file spinor_ops.h
   @brief prototype functions for various spinor operations
 */
#ifndef SPINOR_OPS_H
#define SPINOR_OPS_H

/**
   @fn void adjoint_spinor( struct spinor *adj , const struct spinor S ) 
   @brief conjugate transpose a spinor
 */
void
adjoint_spinor( struct spinor *adj ,
		const struct spinor S ) ;

/**
   @fn double complex bilinear_trace( const struct spinor A , const struct spinor B )
   @brief trace of the product of two spinors Tr[ A*B ]
 */
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B ) ;

/**
   @fn void full_adj( struct spinor *adj , const struct spinor S , const struct gamma G5 )
   @brief computes gamma_5 adj( S ) gamma_5, puts result in adj
 */
void
full_adj( struct spinor *adj ,
	  const struct spinor S ,
	  const struct gamma G5 ) ;

/**
   @fn void gamma_mul_l( struct spinor *res , const struct gamma GAMMA )
   @brief atomic left multiply a spinor by a gamma matrix
 */
void
gamma_mul_l( struct spinor *res ,
	     const struct gamma GAMMA ) ;

/**
   @fn void gamma_mul_r( struct spinor *res , const struct gamma GAMMA )
   @brief atomic right multiply a spinor by a gamma matrix
 */
void
gamma_mul_r( struct spinor *res ,
	     const struct gamma GAMMA ) ;

/**
   @fn void gamma_mul_lr( struct spinor *S , const struct gamma GLEFT , const struct gamma GRIGHT )
   @brief multiply a spinor from the left and right by gamma matrices
 */
void
gamma_mul_lr( struct spinor *S , 
	      const struct gamma GLEFT ,
	      const struct gamma GRIGHT ) ;

/**
   @fn void gauge_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
gauge_spinor( struct spinor *res ,  
	      const double complex link[ NCNC ] ,
	      const struct spinor S ) ;

/**
   @fn void gaugedag_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a daggered link res = link^{\dagger} S
 */
void
gaugedag_spinor( struct spinor *res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

/**
   @fn void gauge_spinordag( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S ) 
   @brief computes U.S^{\dagger}
 */
void
gauge_spinordag( struct spinor *res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S ) ;

/**
   @fn double complex meson_contract( const struct gamma GSNK ,	const struct spinor bwd , const struct gamma GSRC , const struct spinor fwd , const struct gamma G5 )
   @brief contraction of bilinears, uses gamma_5 hermiticity
   @params GSNK :: sink gamma matrix
   @params bwd :: backward propagator solution
   @params GSRC :: source gamma matrix
   @params fwd :: forward propagator solution
   @params G5 :: gamma matrix #GAMMA_5

   Performs \[ \text{Tr}\left[ \gamma_{SNK} \gamma_5 ( \text{bwd} )^{\dagger} \gamma_5 \gamma_{SRC} \text{fwd} \right] \]

   @return the trace
 */
double complex
meson_contract( const struct gamma GSNK ,		
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 ) ;

/**
   @fn void gauge_spinor( struct spinor *res , const double complex link[ NCNC ] , const struct spinor S )
   @brief multiplies a spinor with a link matrix res = link * S
 */
void
spinor_gauge( struct spinor *res ,
	      const struct spinor S ,
	      const double complex link[ NCNC ] ) ;

/**
   @fn void spinordag_gauge( struct spinor *res , const struct spinor S , const double complex link[ NCNC ] )
   @brief computes S^{\dagger}.U
 */
void
spinordag_gauge( struct spinor *res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

/**
   @fn void spinor_gaugedag( struct spinor *res , const struct spinor S , const double complex link[ NCNC ] )
   @brief computes S.U^{\dagger}
 */
void
spinor_gaugedag( struct spinor *res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] ) ;

#endif
