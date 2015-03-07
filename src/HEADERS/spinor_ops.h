/**
   @file spinor_ops.h
   @brief prototype functions for various spinor operations
 */
#ifndef SPINOR_OPS_H
#define SPINOR_OPS_H

/**
   @fn void adjoint_spinor( struct spinor *__restrict adj , const struct spinor S ) 
   @brief conjugate transpose a spinor
 */
void
adjoint_spinor( struct spinor *__restrict adj ,
		const struct spinor S ) ;

/**
   @fn double complex bilinear_trace( const struct spinor A , const struct spinor B )
   @brief trace of the product of two spinors Tr[ A*B ]
 */
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B ) ;

/**
   @fn void chiral_to_nrel( struct spinor *S )
   @brief convert chiral basis to non-relativistic

   change from chiral into nrel basis
    
   chiral \f$ S = T * S * T^\dagger \f$
   here T:
   
           1/s    0       1/s     0
           0      1/s     0       1/s
           -1/s   0       1/s     0
           0      -1/s    0       1/s
    
   with \f$s = \sqrt{2} \f$
 */
void
chiral_to_nrel( struct spinor *S ) ;

/**
   @fn void full_adj( struct spinor *__restrict adj , const struct spinor S , const struct gamma G5 )
   @brief computes \f$ gamma_5 adj( S ) gamma_5 \f$ , puts result in adj
 */
void
full_adj( struct spinor *__restrict adj ,
	  const struct spinor S ,
	  const struct gamma G5 ) ;

/**
   @fn void gamma_mul_l( struct spinor *__restrict res , const struct gamma GAMMA )
   @brief atomic left multiply a spinor by a gamma matrix
 */
void
gamma_mul_l( struct spinor *__restrict res ,
	     const struct gamma GAMMA ) ;

/**
   @fn void gamma_mul_r( struct spinor *__restrict res , const struct gamma GAMMA )
   @brief atomic right multiply a spinor by a gamma matrix
 */
void
gamma_mul_r( struct spinor *__restrict res ,
	     const struct gamma GAMMA ) ;

/**
   @fn void gamma_mul_lr( struct spinor *__restrict S , const struct gamma GLEFT , const struct gamma GRIGHT )
   @brief multiply a spinor from the left and right by gamma matrices
 */
void
gamma_mul_lr( struct spinor *__restrict S , 
	      const struct gamma GLEFT ,
	      const struct gamma GRIGHT ) ;

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
   @fn double complex meson_contract( const struct gamma GSNK ,	const struct spinor bwd , const struct gamma GSRC , const struct spinor fwd , const struct gamma G5 )
   @brief contraction of bilinears, uses gamma_5 hermiticity
   @param GSNK :: sink gamma matrix
   @param bwd :: backward propagator solution
   @param GSRC :: source gamma matrix
   @param fwd :: forward propagator solution
   @param G5 :: gamma 5

   @return \f[
    \textrm{Tr} \left[ \gamma_{\textrm{GSNK}} \gamma_5 ( \textrm{bwd} )^{\dagger} \gamma_5 \gamma_{\textrm{GSRC}} ( \textrm{fwd} ) \right]
   \f]
 */
double complex
meson_contract( const struct gamma GSNK ,		
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 ) ;

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

#endif
