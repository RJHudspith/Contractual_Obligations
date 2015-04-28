/**
   @file matrix_ops.h
   @brief prototype functions for general matrix operations
 */
#ifndef MATRIX_OPS_SSE_H
#define MATRIX_OPS_SSE_H

#ifdef HAVE_EMMINTRIN_H

/**
   @fn void add_mat( __m128d *a , const __m128d *b )
   @brief atomically add color matrix b to a
 */
void
add_mat( __m128d *a ,
	 const __m128d *b ) ;

/**
   @fn inline void colormatrix_equiv( __m128d *a , const __m128d *b )
   @brief equate two colormatrices
 */
void
colormatrix_equiv( double complex *a ,
		   const double complex *b ) ;

/**
   @fn inline void colormatrix_equiv_f2d( double complex a[ NCNC ] , const float complex b[ NCNC ] )
   @brief cast from single to working (double) precision a color matrix
 */
void
colormatrix_equiv_f2d( double complex a[ NCNC ] ,
		       const float complex b[ NCNC ] ) ;


/**
   @fn __m128d colortrace( const __m128d *a )
   @brief take the color trace of a matrix
 */
__m128d
colortrace( const __m128d *a ) ;

/**
   @fn __m128d colortrace_prod( const __m128d *a , const __m128d *b )
   @brief trace of the product of two color matrices
 */
__m128d
colortrace_prod( const __m128d *a , 
		 const __m128d *b ) ;

/**
   @fn constant_mul_gauge( double complex *res , const double complex constant , const double complex *U )
   @brief computes \f$ res = constant * U \f$
 */
void
constant_mul_gauge( double complex *res , 
		    const double complex constant ,
		    const double complex *U ) ;

/**
   @fn void dagger_gauge( __m128d *res , const __m128d *U )
   @brief computes \f$ res = U^{\dagger} \f$
 */
void
dagger_gauge( __m128d *res ,
	      const __m128d *U ) ;

/**
   @fn void multab( __m128d *__restrict a , const __m128d *__restrict b , const __m128d *__restrict c )
   @brief performs the matrix multiplication \f$ a = b * c \f$
 */
void 
multab( __m128d *__restrict a , 
	const __m128d *__restrict b , 
	const __m128d *__restrict c ) ;

/**
   @fn void multabdag( __m128d *__restrict a , const __m128d *__restrict b , const __m128d *__restrict c )
   @brief performs the matrix multiplication \f$ a = b^{\dagger} * c \f$
 */
void 
multabdag( __m128d *__restrict a , 
	   const __m128d *__restrict b , 
	   const __m128d *__restrict c ) ;

/**
   @fn void multab_dag( __m128d *__restrict a , const __m128d *__restrict b , const __m128d *__restrict c )
   @brief performs the matrix multiplication \f$ a = b * c ^{\dagger} \f$
 */
void 
multab_dag( __m128d *__restrict a , 
	    const __m128d *__restrict b , 
	    const __m128d *__restrict c ) ;

/**
   @fn void multab_dagdag( __m128d *__restrict a , const __m128d *__restrict b , __m128d *__restrict c )
   @brief performs the matrix multiplication \f$ a = b^{\dagger} * c ^{\dagger} \f$
 */
void 
multab_dagdag( __m128d *__restrict a , 
	       const __m128d *__restrict b , 
	       const __m128d *__restrict c ) ;

/**
   @fn void print_colormatrix( const double complex a[ NCNC ] )
   @brief print to stdout the color matrix a
 */
void
print_colormatrix( const double complex a[ NCNC ] ) ;

#endif

#endif
