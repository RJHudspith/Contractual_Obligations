/**
   @file bar_projections.h
   @brief baryon projection definitions
 */
#ifndef BAR_PROJECTIONS_H
#define BAR_PROJECTIONS_H

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

/**
   @fn void compute_p_psq( double p[ ND ] , double *p2 , const struct veclist momentum )
   @brief compute momenta and psq
 */
void
compute_p_psq( double p[ ND ] ,
	       double *p2 ,
	       const struct veclist momentum ) ;

/**
   @fn void P32( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin 3/2 projection
 **/
void
P32( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn void P11( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin 1/2 (11) projection
 **/
void
P11( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn void P12( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin 1/2 (12) projection
 **/
void
P12( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn void P21( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin 1/2 (21) projection
 **/
void
P21( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn void P22( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin 1/2 (22) projection
 **/
void
P22( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn void P00( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief null matrix projection
 **/
void
P00( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx ) ;

/**
   @fn int spinproj( double complex *Gik , void (*p)( double complex *proj , const size_t i , const size_t j , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx ) , const struct mcorr **corr , const size_t i , const size_t k , const size_t t , const struct gamma *GAMMA , const struct veclist *momentum , const size_t pidx )
   @brief spin project a baryon correlation function
 */
int
spinproj( double complex *Gik ,
	  void (*p)( double complex *proj ,
		     const size_t i ,
		     const size_t j ,
		     const struct gamma *GAMMA ,
		     const struct veclist *momentum ,
		     const size_t pidx ) ,
	  const struct mcorr **corr ,
	  const size_t i ,
	  const size_t k ,
	  const size_t t ,
	  const struct gamma *GAMMA ,
	  const struct veclist *momentum ,
	  const size_t pidx ) ;

#endif
