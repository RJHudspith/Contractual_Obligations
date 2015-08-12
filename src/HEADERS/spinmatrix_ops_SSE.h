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
   @fn double complex gammaspinmatrix_trace( const struct gamma G , const void *spinmatrix )
   @brief performs the spin trace of \gamma * spinmatrix
 */
double complex
gammaspinmatrix_trace( const struct gamma G ,
		       const void *spinmatrix ) ;

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
   @fn void zero_spinmatrix( void *spinmatrix )
   @brief zero a spinmatrix
 */
void
zero_spinmatrix( void *spinmatrix ) ;

#endif

#endif
