/**
   @file contractions.h
   @brief traces, brute force gamma multiplies and adjoints
 */

#ifndef CONTRACTIONS_SSE_H
#define CONTRACTIONS_SSE_H

#ifdef HAVE_EMMINTRIN_H

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
   @fn double complex simple_meson_contract( const struct gamma GSNK , const struct spinor bwd , const struct gamma GSRC , const struct spinor fwd )
   @brief does a simple meson contraction
   
   @return \f[
   \textrm{Tr} \left[ \gamma_{\textrm{GSNK}} \textrm{bwd} \gamma_{\textrm{GSRC}} ( \textrm{fwd} ) \right]
   \f]
 */
double complex
simple_meson_contract( const struct gamma GSNK ,		
		       const struct spinor bwd , 
		       const struct gamma GSRC ,
		       const struct spinor fwd ) ;

#endif

#endif
