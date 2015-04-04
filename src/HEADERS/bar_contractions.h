/**
   @file bar_contractions.h
   @brief prototype functions for baryon contractions
 */

#ifndef BAR_CONTRACTIONS_H
#define BAR_CONTRACTIONS_H

/**
   @fn void baryon_contract_site( )
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
   @fn void baryon_contract_omega_site( )
   @brief performs the baryon contraction at a site and has the 

 */
void
baryon_contract_site( double complex **term ,
		      const struct spinor S1 , 
		      const struct spinor S2 , 
		      const struct spinor S3 , 
		      const struct gamma Cgmu ,
		      const struct gamma CgmuT ) ;

#endif
