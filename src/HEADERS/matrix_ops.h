/**
   @file matrix_ops.h
   @brief prototype functions for general matrix operations
 */
#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

/**
   @fn inline void add_mat( double complex *a , const double complex *b )
   @brief atomically add color matrix b to a
 */
inline void
add_mat( double complex *a ,
	 const double complex *b ) ;

/**
   @fn inline void colormatrix_equiv( double complex a[ NCNC ] , const double complex b[ NCNC ] )
   @brief equate two colormatrices
 */
inline void
colormatrix_equiv( double complex a[ NCNC ] ,
		   const double complex b[ NCNC ] ) ;

/**
   @fn inline void colormatrix_equiv_f2d( double complex a[ NCNC ] , const float complex b[ NCNC ] )
   @brief cast from single to working (double) precision a color matrix
 */
inline void
colormatrix_equiv_f2d( double complex a[ NCNC ] ,
		       const float complex b[ NCNC ] ) ;

/**
   @fn inline double complex colortrace_prod( const double complex *a , const double complex *b )
   @brief trace of the product of two color matrices
 */
inline double complex
colortrace_prod( const double complex *a , 
		 const double complex *b ) ;

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
