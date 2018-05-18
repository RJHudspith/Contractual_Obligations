/**
   @file mmul.h
   @brief prototype declarations for color matrix multiplies
 */
#ifndef MMUL_H
#define MMUL_H

#ifdef HAVE_EMMINTRIN_H

#include "mmul_SSE.h"

#else

/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b*c
 **/
void 
multab( double complex a[ NCNC ] , 
	const double complex b[ NCNC ] , 
	const double complex c[ NCNC ] ) ;

/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b*c, b and c are members of SU(#NC)
 **/
void 
multab_suNC( double complex a[ NCNC ] , 
	     const double complex b[ NCNC ] , 
	     const double complex c[ NCNC ] ) ;
/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b^\dagger * c
 **/
void 
multabdag( double complex a[ NCNC ] , 
	   const double complex b[ NCNC ] , 
	   const double complex c[ NCNC ] ) ;
/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b^\dagger * c b and c are members of SU(#NC)
 **/
void 
multabdag_suNC( double complex a[ NCNC ] , 
		const double complex b[ NCNC ] , 
		const double complex c[ NCNC ] ) ;
/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b*c^\dagger
 **/
void 
multab_dag( double complex a[ NCNC ] , 
	    const double complex b[ NCNC ] , 
	    const double complex c[ NCNC ] ) ;
/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b*c^\dagger b and c are members of SU(#NC)
 **/
void 
multab_dag_suNC( double complex a[ NCNC ] , 
		 const double complex b[ NCNC ] , 
		 const double complex c[ NCNC ] ) ;
/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b^\dagger * c^\dagger
 **/
void 
multab_dagdag( double complex a[ NCNC ] , 
	       const double complex b[ NCNC ] , 
	       const double complex c[ NCNC ] ) ;

/**
   @fn void multab( double complex a[ NCNC ] , const double complex b[ NCNC ] , const double complex c[ NCNC ] ) ;
   @brief #NC x #NC matrix multiply a = b^dagger *c^dagger b and c are members of SU(#NC)
 **/
void 
multab_dagdag_suNC( double complex a[ NCNC ] , 
		    const double complex b[ NCNC ] , 
		    const double complex c[ NCNC ] ) ;

#endif

#endif
