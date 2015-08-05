/**
   @file baryon_projections.h
   @brief prototype declarations for baryon projections
 */
#ifndef BARYON_PROJECTIONS_H
#define BARYON_PROJECTIONS_H

// spin 1/2 projection types
typedef enum { NONE , OneHalf_11 , OneHalf_12 , OneHalf_21 ,
	       OneHalf_22 , ThreeHalf } spinhalf ;


// enum for the projections we allow
typedef enum { L0 , L1 , L2 , L3 , L4 , L5 } bprojection ;

/**
   @fn double complex* baryon_project( const struct mcorr **corr , const struct gamma *GAMMA , const struct veclist *momentum , const size_t GSRC , const size_t GSNK , const size_t p , const bprojection parity_proj , const spinhalf spin_proj ) 
   @brief project out the open dirac indices
 */
double complex*
baryon_project( const struct mcorr **corr ,
		const struct gamma *GAMMA ,
		const struct veclist *momentum ,
		const size_t GSRC ,
		const size_t GSNK ,
		const size_t p ,
		const bprojection parity_proj ,
		const spinhalf spin_proj ) ;

#endif
