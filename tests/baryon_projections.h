/**
   @file baryon_projections.h
   @brief prototype declarations for baryon projections
 */
#ifndef BARYON_PROJECTIONS_H
#define BARYON_PROJECTIONS_H

// spin 1/2 projection types
typedef enum { OneHalf_11 , OneHalf_12 , OneHalf_21 ,
	       OneHalf_22 , ThreeHalf } spinhalf ;


// enum for the projections we allow
typedef enum { LO , L1 , L2 , L3 , L4 , L5 } bprojection ;

/**
   @fn void spin_project( double complex *Oi , const struct mcorr **corr , const struct gamma *GAMMA , const struct veclist *momentum , const size_t GSRC , const size_t GSNK , const size_t p , const size_t t , const spinhalf projection ) 
   @brief spin-1/2 and 3/2 projections
 */
void
spin_project( double complex *Oi , 
	      const struct mcorr **corr ,
	      const struct gamma *GAMMA ,
	      const struct veclist *momentum ,
	      const size_t GSRC ,
	      const size_t GSNK ,
	      const size_t p ,
	      const size_t t ,
	      const spinhalf projection ) ;

/**
   @fn double complex* baryon_project( const struct mcorr **corr , const struct gamma *GAMMA , const struct veclist *momentum , const size_t GSRC , const size_t GSNK , const size_t p , const bprojection projection ) 
   @brief project out the open dirac indices
 */
double complex*
baryon_project( const struct mcorr **corr ,
		const struct gamma *GAMMA ,
		const struct veclist *momentum ,
		const size_t GSRC ,
		const size_t GSNK ,
		const size_t p ,
		const bprojection projection ) ;

#endif
