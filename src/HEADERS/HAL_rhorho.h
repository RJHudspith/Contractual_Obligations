/**
   @file HAL_rhorho.h
   @brief prototype declarations for HAL-QCD approach to rho-rho DM
 */
#ifndef HAL_RHORHO_H
#define HAL_RHORHO_H

/**
   @fn int HALrhorho_contract( double complex *in , double complex *out , fftw_plan forward , fftw_plan backward , struct spinmatrix **blk11 , struct spinmatrix **blk12 , struct spinmatrix **blk21 , struct spinmatrix **blk22 , const struct spinor *S , const struct gamma *GAMMAS , const size_t GSRC , const size_t GSNK )
   @brief do the HAL thing for baryons
 */
int
HALrhorho_contract( double complex *in ,
		    double complex *out ,
		    fftw_plan forward ,
		    fftw_plan backward ,
		    struct spinmatrix **blk11 ,
		    struct spinmatrix **blk12 ,
		    struct spinmatrix **blk21 ,
		    struct spinmatrix **blk22 ,
		    const struct spinor *S ,
		    const struct gamma *GAMMAS ,
		    const size_t GSRC ,
		    const size_t GSNK ,
		    const int nmom[1] ,
		    const struct veclist *list ) ;


#endif
