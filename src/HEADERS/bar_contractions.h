/**
   @file bar_contractions.h
   @brief prototype functions for baryon contractions
 */

#ifndef BAR_CONTRACTIONS_H
#define BAR_CONTRACTIONS_H

/**
   @fn void baryon_contract_site( double complex **term , const struct spinor S1 , const struct spinor S2 , const struct spinor S3 , const struct gamma Cgmu , const struct gamma CgmuT )
   @brief performs the baryon contraction at a site and has the 
 */
void
baryon_contract_site( double complex **term ,
		      const struct spinor S1 , 
		      const struct spinor S2 , 
		      const struct spinor S3 , 
		      const struct gamma Cgmu ,
		      const struct gamma CgmuT ) ;

/**
   @fn void baryon_contract_omega_site( double complex **term , const struct spinor S1 , const struct spinor S2 , const struct spinor S3 , const struct gamma Cgmu , const struct gamma CgmuT )
   @brief performs the baryon contraction at a site and has the 
 */
void
baryon_contract_omegasite( double complex **term ,
			   const struct spinor S1 , 
			   const struct spinor S2 , 
			   const struct spinor S3 , 
			   const struct gamma Cgmu ,
			   const struct gamma CgmuT ) ;

/**
   @fn void baryon_contract_site_mom( double complex **in , const struct spinor S1 , const struct spinor S2 , const struct spinor S3 , const struct gamma Cgmu , const struct gamma CgmuT , const int GSRC , const int site )
   @brief perform a baryon contraction accumulating site-wise result into flat "in" array for FFT-ing
 */
void
baryon_contract_site_mom( double complex **in ,
			  const struct spinor S1 , 
			  const struct spinor S2 , 
			  const struct spinor S3 , 
			  const struct gamma Cgmu ,
			  const struct gamma CgmuT ,
			  const int GSRC ,
			  const int site ) ;

#endif
