/**
   @file bar_contractions.h
   @brief prototype functions for baryon contractions
 */
#ifndef BAR_CONTRACTIONS_H
#define BAR_CONTRACTIONS_H

/**
   @fn void baryon_contract_omega_site( double complex **term , const struct spinor S1 , const struct spinor S2 , const struct spinor S3 , const struct gamma Cgmu , const struct gamma GgmuD )
   @brief performs the baryon contraction at a site and has the 
 */
void
baryon_contract_omega_site( double complex **term ,
			    const struct spinor S1 , 
			    const struct spinor S2 , 
			    const struct spinor S3 , 
			    const struct gamma Cgmu ,
			    const struct gamma GgmuD ) ;

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
   @fn void baryon_contract_site_mom( double complex **in , const struct spinor S1 , const struct spinor S2 , const struct spinor S3 , const struct gamma Cgmu , const struct gamma CgmuT , const size_t GSRC , const size_t site )
   @brief perform a baryon contraction accumulating site-wise result into flat "in" array for FFT-ing
 */
void
baryon_contract_site_mom( double complex **in ,
			  const struct spinor S1 , 
			  const struct spinor S2 , 
			  const struct spinor S3 , 
			  const struct gamma Cgmu ,
			  const struct gamma CgmuT ,
			  const size_t GSRC ,
			  const size_t site ) ;

/**
   @fn void baryon_contract_walls( struct mcorr **corr , const struct spinor SUM1 , const struct spinor SUM2 , const struct spinor SUM3 , const struct gamma *GAMMAS , const size_t t , const baryon_type btype )
   @brief perform the Wall-Wall Baryon contractions
   @warning must be called inside the parallel environment
 */
void
baryon_contract_walls( struct mcorr **corr , 
		       const struct spinor SUM1 ,
		       const struct spinor SUM2 ,
		       const struct spinor SUM3 ,
		       const struct gamma *GAMMAS ,
		       const size_t t ,
		       const baryon_type btype ) ;

/**
   @fn void baryon_momentum_project( struct measurements *M , const size_t stride1 , const size_t stride2 , const size_t t , const baryon_type btype )
   @brief perform the momentum projection for our baryons
 */
void
baryon_momentum_project( struct measurements *M ,
			 const size_t stride1 , 
			 const size_t stride2 ,
			 const size_t t ,
			 const baryon_type btype ) ;

#endif
