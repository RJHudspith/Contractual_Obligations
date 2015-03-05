/**
   @file matrix_ops.h
   @brief prototype functions for general matrix operations
 */
#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

/**
   @fn inline double complex colortrace_prod( double complex *a , double complex *b )
   @brief trace of the product of two color matrices
 */
inline double complex
colortrace_prod( double complex *a , 
		 double complex *b ) ;

/**
   @fn constant_mul_gauge( double complex *res , const double complex constant , const double complex *U )
   @brief computes \f$ res = constant * U \f$
 */
inline void
constant_mul_gauge( double complex *res , 
		    const double complex constant ,
		    const double complex *U ) ;

/**
   @fn inline void dagger_gauge( double complex *res , const double complex *U )
   @brief computes \f$ res = U^{\dagger} \f$
 */
inline void
dagger_gauge( double complex *res ,
	      const double complex *U ) ;

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

#endif
