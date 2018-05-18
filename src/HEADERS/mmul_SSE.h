#ifndef MMUL_SSE_H
#define MMUL_SSE_H

void 
multab( __m128d *__restrict a , 
	const __m128d *__restrict b , 
	const __m128d *__restrict c ) ;

void 
multab_suNC( __m128d *__restrict a , 
	     const __m128d *__restrict b , 
	     const __m128d *__restrict c ) ;

void 
multabdag( __m128d *__restrict a , 
	   const __m128d *__restrict b , 
	   const __m128d *__restrict c ) ;

void 
multabdag_suNC( __m128d *__restrict A , 
		const __m128d *__restrict B , 
		const __m128d *__restrict C ) ;

void 
multab_dag_suNC( __m128d *__restrict A , 
		 const __m128d *__restrict B , 
		 const __m128d *__restrict C ) ;

void 
multab_dag( __m128d *__restrict a , 
	    const __m128d *__restrict b , 
	    const __m128d *__restrict c ) ;

void 
multab_dag_suNC( __m128d *__restrict a , 
		 const __m128d *__restrict b , 
		 const __m128d *__restrict c ) ;

void 
multab_dagdag( __m128d *__restrict a , 
	       const __m128d *__restrict b , 
	       const __m128d *__restrict c ) ;

void 
multab_dagdag_suNC( __m128d *__restrict A , 
		    const __m128d *__restrict B , 
		    const __m128d *__restrict C ) ;
#endif
