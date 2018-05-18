/**
   @file matrix_ops.h
   @brief prototype functions for general matrix operations
 */
#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#ifdef HAVE_EMMINTRIN_H

#include "matrix_ops_SSE.h"

#else

/**
   @fn void add_mat( double complex *__restrict a , const double complex *__restrict b )
   @brief atomically add color matrix b to a
 */
void
add_mat( double complex *__restrict a ,
	 const double complex *__restrict b ) ;

/**
   @fn void colormatrix_equiv( double complex *__restrict a , const double complex *__restrict b )
   @brief equate two colormatrices
 */
void
colormatrix_equiv( double complex *__restrict a ,
		   const double complex *__restrict b ) ;

/**
   @fn void colormatrix_equiv_f2d( double complex a[ NCNC ] , const float complex b[ NCNC ] )
   @brief cast from single to working (double) precision a color matrix
 */
void
colormatrix_equiv_f2d( double complex a[ NCNC ] ,
		       const float complex b[ NCNC ] ) ;

/**
   @fn void colormatrix_equiv_d2f( float complex a[ NCNC ] , const double complex b[ NCNC ] )
   @brief cast from working (double) precision to a single precision color matrix
 */
void
colormatrix_equiv_d2f( float complex a[ NCNC ] ,
		       const double complex b[ NCNC ] ) ;

/**
   @fn void colormatrix_iSaxpy( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double S )
   @brief computes a += S*(i*b) where a and b are colormatrices and S is a scalar
 */
void
colormatrix_iSaxpy( double complex a[ NCNC ] ,
		    const double complex b[ NCNC ] ,
		    const double S ) ;

/**
   @fn void colormatrix_Saxpy( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double S )
   @brief computes a += S*b where a and b are colormatrices and S is a scalar
 */
void
colormatrix_Saxpy( double complex a[ NCNC ] ,
		   const double complex b[ NCNC ] ,
		   const double S ) ;

/**
   @fn void colormatrix_Sa_xmy( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double S )
   @brief computes a[i] = S*( a[i] - b[i] )
 */
void
colormatrix_Sa_xmy( double complex a[ NCNC ] ,
		    const double complex b[ NCNC ] ,
		    const double S ) ;

/**
   @fn double complex colortrace( const double complex a[ NCNC ] )
   @brief compute the trace of a color matrix
 */
double complex
colortrace( const double complex a[ NCNC ] ) ;

/**
   @fn double complex colortrace_prod( const double complex *__restrict a , const double complex *__restrict b )
   @brief trace of the product of two color matrices
 */
double complex
colortrace_prod( const double complex *__restrict a , 
		 const double complex *__restrict b ) ;

/**
   @fn void constant_mul_gauge( double complex *__restrict res , const double complex constant , const double complex *__restrict U )
   @brief computes \f$ res = constant * U \f$
 */
void
constant_mul_gauge( double complex *__restrict res , 
		    const double complex constant ,
		    const double complex *__restrict U ) ;

/**
   @fn void dagger_gauge( double complex *__restrict res , const double complex *__restrict U )
   @brief computes \f$ res = U^{\dagger} \f$
 */
void
dagger_gauge( double complex *__restrict res ,
	      const double complex *__restrict U ) ;

/**
   @fn double complex LU_det( const int N , const double complex U[ N*N ] )
   @brief determinant from an LU factorisation
   @return the determinant of U
 */
double complex
LU_det( const int N , 
	const double complex U[ N*N ] ) ;

/**
   @fn void print_colormatrix( const double complex a[ NCNC ] )
   @brief print to stdout the color matrix a
 */
void
print_colormatrix( const double complex a[ NCNC ] ) ;

/**
   @fn void zero_colormatrix( const double complex a[ NCNC ] )
   @brief zero a colormatrix
 */
void
zero_colormatrix( double complex a[ NCNC ] ) ;

#endif

#endif
