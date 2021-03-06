/**
   @file spinmatrix_ops_SSE.h
   @brief SSEd spinmatrix linear algebra
 */
#ifndef SPINMATRIX_OPS_SSE_H
#define SPINMATRIX_OPS_SSE_H

#ifdef HAVE_EMMINTRIN_H

/**
   @fn void atomic_add_spinmatrices( void *res , const void *D )
   @brief atomically spinmatrices res += D
 */
void
atomic_add_spinmatrices( void *res ,
			 const void *D ) ;

/**
   @fn void compute_pslash( void *pslash , const struct gamma *GAMMA , const double p[ ND ] )
   @brief compute pslash \f$ \sum_{i=0}^{i<DIMS} p_i.\gamma_i \f$
 */
void
compute_pslash( void *pslash , 
		const struct gamma *GAMMA ,
		const double p[ ND ] ) ;

/**
   @fn void gamma_spinmatrix( void *spinmatrix , const struct gamma G ) 
   @brief atomically left multiply by a gamma, spinmatrix = G * spinmatrix
 */
void
gamma_spinmatrix( void *spinmatrix ,
		  const struct gamma G ) ;

/**
   @fn void gamma_spinmatrix_lr( void *spinmatrix , const struct gamma GLEFT , const struct gamma GRIGHT ) 
   @brief atomically left and right multiply by a gamma, spinmatrix = GLEFT * spinmatrix * GRIGHT
 **/
void
gamma_spinmatrix_lr( struct spinmatrix *S ,
		     const struct gamma GLEFT ,
		     const struct gamma GRIGHT ) ;

/**
   @fn double complex gammaspinmatrix_trace( const struct gamma G , const void *spinmatrix )
   @brief performs the spin trace of \gamma * spinmatrix
 */
double complex
gammaspinmatrix_trace( const struct gamma G ,
		       const void *spinmatrix ) ;

/**
   @fn void get_spinmatrix( void *spinmatrix , const struct spinor S , const size_t c1 , const size_t c2 )
   @brief from our spinor struct get a matrix with only spin indices
 */
void
get_spinmatrix( void *spinmatrix , 
		const struct spinor S ,
		const size_t c1 ,
		const size_t c2 ) ;

/**
   @fn void identity_spinmatrix( void *spinmatrix )
   @brief set spinmatrix to the identity
 */
void
identity_spinmatrix( void *spinmatrix ) ;

/**
   @fn void spinmatrix_gamma( void *spinmatrix , const struct gamma G ) 
   @brief right multiply by a gamma, res = spinmatrix * G
 */
void
spinmatrix_gamma( void *spinmatrix ,
		  const struct gamma G ) ;

/**
   @fn void spinmatrix_multiply( void *a , const void *b , const void *c )
   @brief multiply two spinmatrices a = b * c
 */
void
spinmatrix_multiply( void *a ,
		     const void *b ,
		     const void *c ) ;

/**
   @fn void spinmatrix_multiply_T( void *a , const void *b , const void *c )
   @brief multiply two spinmatrices where the second is transposed a = b * c^T
*/
void
spinmatrix_multiply_T( void *a ,
		       const void *b ,
		       const void *c ) ;

/**
   @fn void spinmatrix_mulconst( void *spinmatrix , const double factor )
   @brief atomically multiply a spinmatrix by a scalar
 */
void
spinmatrix_mulconst( void *spinmatrix , 
		     const double factor ) ;


/**
   @fn double complex spinmatrix_trace( const void *spinmatrix )
   @brief trace of a spinmatrix
 */
double complex
spinmatrix_trace( const void *spinmatrix ) ;

/**
   @fn double complex trace_prod_spinmatrices( const void *a , const void *b )
   @brief trace of the product of two spinmatrices
   @return Tr[a.b]
 */
double complex
trace_prod_spinmatrices( const void *a , 
			 const void *b ) ;

/**
   @fn double complex trace_prod_spinmatrices_dag( const void *a , const void *b )
   @brief trace of the product of two spinmatrices
   @return Tr[a.b^\dag]
 */
double complex
trace_prod_spinmatrices_dag( const void *a , 
			     const void *b ) ;

/**
   @fn transpose_spinmatrix( void *a )
   @brief transpose a spinmatrix
*/
void
transpose_spinmatrix( void *a ) ;

/**
   @fn void zero_spinmatrix( void *spinmatrix )
   @brief zero a spinmatrix
 */
void
zero_spinmatrix( void *spinmatrix ) ;

#endif

#endif
