/**
   @file brutal_mesons.h
   @brief brute-force meson contraction code for sanity testing
 */

#ifndef BRUTAL_MESONS_H
#define BRUTAL_MESONS_H

/**
   @fn double complex meson_trace( const struct gamma GSNK , const struct spinor S2 , const struct gamma GSRC , const struct spinor S1 ) 
   @brief brute force meson trace
 */
double complex
meson_trace( const struct gamma GSNK ,
	     const struct spinor S2 ,
	     const struct gamma GSRC ,
	     const struct spinor S1 ) ;

/**
   @fn int single_mesons_bruteforce( FILE *prop1 , const proptype proptype1 , const char *outfile )
   @brief same as meson contraction code, only more brutal
 */
int
single_mesons_bruteforce( FILE *prop1 , 
			  const proptype proptype1 ,
			  const char *outfile ) ;

#endif
