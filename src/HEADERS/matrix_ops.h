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
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication \f$ a = b * c \f$
 */
void 
multab( double complex a[ NCNC ] , 
	const double complex b[ NCNC ] , 
	const double complex c[ NCNC ] ) ;

/**
   @fn void multabdag( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication \f$ a = b^{\dagger} * c \f$
 */
void 
multabdag( double complex a[ NCNC ] , 
	   const double complex b[ NCNC ] , 
	   const double complex c[ NCNC ] ) ;

/**
   @fn void multab_dag( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication \f$ a = b * c ^{\dagger} \f$
 */
void 
multab_dag( double complex a[ NCNC ] , 
	    const double complex b[ NCNC ] , 
	    const double complex c[ NCNC ] ) ;

/**
   @fn void multab_dagdag( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] )
   @brief performs the matrix multiplication \f$ a = b^{\dagger} * c ^{\dagger} \f$
 */
void 
multab_dagdag( double complex a[ NCNC ] , 
	       const double complex b[ NCNC ] , 
	       const double complex c[ NCNC ] ) ;

/**
   @fn void print_colormatrix( const double complex a[ NCNC ] )
   @brief print to stdout the color matrix a
 */
void
print_colormatrix( const double complex a[ NCNC ] ) ;

#endif

#endif
