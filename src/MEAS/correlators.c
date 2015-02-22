/**
   @file correlators.c
   @brief correlation function calculations

   TODO :: Profile this -> see if we can gain much by unrolling the c1, c2 loop
   could be done with double pointers and some shenanigans I guess ..
 */

#include "common.h"

#include "gammas.h"

// singlet has to be real so no point in computing imaginary part
double complex 
local_meson_correlator_singlet( const struct spinor S1 , 
				const struct gamma ADJ ,
				const struct gamma SRC ,
				const struct gamma SNK )
{
  // loop counters
  int d1 , d2 , c1 , c2 ;
  int sign ; // permutation of the fourth roots of unity 

  register double complex corr = 0.0 ;
  register double rloc_corr , iloc_corr ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {

      // map the gamma multiplied indices for legibility 
      const int id1 = SRC.ig[ ADJ.ig[ d1 ] ] ;
      const int id2 = SNK.ig[ ADJ.ig[ d2 ] ] ;

      // adjust the sign is a permutation of the fourth roots of unity for all sensible
      // gamma matrix conventions
      // only the 1 and -1 contribute so we just keep those
      sign = ( ( SRC.g[d1] + ADJ.g[id1] + SNK.g[d2] + ADJ.g[id2] ) & 3 ) ;
     
      rloc_corr = iloc_corr = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  //loc_corr += conj( S1.D[d1][d2].C[c1][c2] ) * S2.D[id1][id2].C[c1][c2];
	  rloc_corr += creal( S1.D[d1][d2].C[c1][c2] ) * creal( S1.D[id1][id2].C[c1][c2] ) + \
	    cimag( S1.D[d1][d2].C[c1][c2] ) * cimag( S1.D[id1][id2].C[c1][c2] ) ;
	  //
	  iloc_corr += creal( S1.D[d1][d2].C[c1][c2] ) * cimag( S1.D[id1][id2].C[c1][c2] ) - \
	    cimag( S1.D[d1][d2].C[c1][c2] ) * creal( S1.D[id1][id2].C[c1][c2] ) ;
	}
      }
      // is just a permutation with this basis this requires NS = 4
      switch( sign ) {
      case 0 : corr += rloc_corr + I * iloc_corr ; break ;
      case 1 : corr += -iloc_corr + I * rloc_corr ; break ;
      case 2 : corr += -rloc_corr - I * iloc_corr ; break ;
      case 3 : corr += iloc_corr - I * rloc_corr ; break ;
      }
    } 
  }
  return corr ;
}

// compute the local meson correlator
double complex 
local_meson_correlator( const struct spinor S1 , 
			const struct spinor S2 , 
			const struct gamma ADJ ,
			const struct gamma SRC ,
			const struct gamma SNK )
{
  // loop counters
  int d1 , d2 , c1 , c2 ;
  int sign ; // permutation of the fourth roots of unity

  register double complex corr = 0.0 ;
  register double rloc_corr , iloc_corr ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {

      // map the gamma multiplied indices for legibility 
      const int id1 = SRC.ig[ ADJ.ig[ d1 ] ] ;
      const int id2 = SNK.ig[ ADJ.ig[ d2 ] ] ;

      // adjust the sign is a permutation of the fourth roots of unity for all sensible
      // gamma matrix conventions
      sign = ( SRC.g[d1] + ADJ.g[id1] + SNK.g[d2] + ADJ.g[id2] )&3 ;

      rloc_corr = 0.0 ;
      iloc_corr = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  //loc_corr += conj( S1.D[d1][d2].C[c1][c2] ) * S2.D[id1][id2].C[c1][c2];
	  rloc_corr += creal( S1.D[d1][d2].C[c1][c2] ) * creal( S2.D[id1][id2].C[c1][c2] ) + \
	    cimag( S1.D[d1][d2].C[c1][c2] ) * cimag( S2.D[id1][id2].C[c1][c2] ) ;
	  //
	  iloc_corr += creal( S1.D[d1][d2].C[c1][c2] ) * cimag( S2.D[id1][id2].C[c1][c2] ) - \
	    cimag( S1.D[d1][d2].C[c1][c2] ) * creal( S2.D[id1][id2].C[c1][c2] ) ;
	}
      }
      // is just a permutation with this basis this requires NS = 4
      switch( sign ) {
      case 0 : corr += rloc_corr + I * iloc_corr ; break ;
      case 1 : corr += -iloc_corr + I * rloc_corr ; break ;
      case 2 : corr += -rloc_corr - I * iloc_corr ; break ;
      case 3 : corr += iloc_corr - I * rloc_corr ; break ;
      }
      //
    } 
  }
  return corr ;
}

// placeholder for conserved-local optimisation
#if 0
if( sign == 2 || sign == 0 ) {
  register double rloc_corr = 0.0 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
#if NC == 3
    rloc_corr += creal( S1.D[d1][d2].C[c1][0] ) * creal( S1.D[id1][id2].C[c1][0] ) ;
    rloc_corr += cimag( S1.D[d1][d2].C[c1][0] ) * cimag( S1.D[id1][id2].C[c1][0] ) ;
    rloc_corr += creal( S1.D[d1][d2].C[c1][1] ) * creal( S1.D[id1][id2].C[c1][1] ) ;
    rloc_corr += cimag( S1.D[d1][d2].C[c1][1] ) * cimag( S1.D[id1][id2].C[c1][1] ) ;
    rloc_corr += creal( S1.D[d1][d2].C[c1][2] ) * creal( S1.D[id1][id2].C[c1][2] ) ;
    rloc_corr += cimag( S1.D[d1][d2].C[c1][2] ) * cimag( S1.D[id1][id2].C[c1][2] ) ;
#else
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      rloc_corr += creal( S1.D[d1][d2].C[c1][c2] ) * creal( S1.D[id1][id2].C[c1][c2] ) + \
	cimag( S1.D[d1][d2].C[c1][c2] ) * cimag( S1.D[id1][id2].C[c1][c2] ) ;
    }
#endif
  }
  corr += ( 1 - sign ) * rloc_corr ;
  //
 }
#endif
