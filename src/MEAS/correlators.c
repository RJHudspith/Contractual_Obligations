/**
   @file correlators.c
   @brief correlation function calculations

   TODO :: Profile this -> see if we can gain much by unrolling the c1, c2 loop
   could be done with double pointers and some shenanigans I guess ..
 */

#include "common.h"

#include "gammas.h"

// compute the local meson correlator
double complex 
local_meson_correlator( const struct spinor S1 , 
			const struct spinor S2 , 
			const int *sadj ,
			const int *iadj, 
			const int mu1 , 
			const int mu2 )
{
  register double complex corr = 0.0 ;
  int imu1[ NS ], imu2[ NS ];
  int sgn1[ NS ] , sgn2[ NS ] ; 
  int id1, id2 , d1 , d2 , c1 , c2 ;
  int sign ; // permutation of the fourth roots of unity

  // Source gamma index 
  gamma_matrix( sgn1, imu1 , mu1 ) ;
  // Sink gamma index 
  gamma_matrix( sgn2 , imu2 , mu2 ) ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {

      // map the gamma multiplied indices for legibility 
      id1 = imu1[iadj[d1]];
      id2 = imu2[iadj[d2]];

      // adjust the sign is a permutation of the fourth roots of unity for all sensible
      // gamma matrix conventions
      sign = ( sgn1[d1] + sadj[id1] + sgn2[d2] + sadj[id2] )&3 ;

      register double rloc_corr = 0.0 , iloc_corr = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  //loc_corr += conj( S1.D[d1][d2].C[c1][c2] ) * S1.D[id1][id2].C[c1][c2];
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
      //
    } 
  }

  return corr ;
}
// this will want the gauge field
double complex 
local_conserved_meson_correlator( const struct spinor S1 , 
				  const struct spinor S2 , 
				  const int *sadj ,
				  const int *iadj , 
				  const int mu1 ,  
				  const int mu2 )
{
  register double complex corr = 0.0 ;
  int sgn1[ NS ], sgn2[ NS ] , sign ;
  int imu1[ NS ] , imu2[ NS ] ;
  int id1, id2 , d1, d2, c1, c2 ;

  // Source gamma index 
  gamma_matrix( sgn1 , imu1 , mu1 ) ;
  // Sink gamma index 
  gamma_matrix( sgn2 , imu2 , mu2 ) ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {

      // map the gamma multiplied indices for legibility 
      id1 = imu1[iadj[d1]];
      id2 = imu2[iadj[d2]];

      // adjust the sign 
      sign = ( sgn1[d1] + sadj[id1] + sgn2[d2] + sadj[id2] )&3 ;

     register double rloc_corr = 0.0 , iloc_corr = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  //loc_corr += conj( S1.D[d1][d2].C[c1][c2] ) * S1.D[id1][id2].C[c1][c2];
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
  return corr;
}
