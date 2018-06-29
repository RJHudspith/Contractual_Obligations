/**
   @file mmul_SSE.h
   @brief vectorised colormatrix multiplies
 */
#ifndef MMUL_SSE_H
#define MMUL_SSE_H

void 
multab( __m128d *__restrict a , 
	const __m128d *__restrict b , 
	const __m128d *__restrict c ) ;

void 
multabdag( __m128d *__restrict a , 
	   const __m128d *__restrict b , 
	   const __m128d *__restrict c ) ;

void 
multab_dag( __m128d *__restrict a , 
	    const __m128d *__restrict b , 
	    const __m128d *__restrict c ) ;

void 
multab_dagdag( __m128d *__restrict a , 
	       const __m128d *__restrict b , 
	       const __m128d *__restrict c ) ;

#endif
