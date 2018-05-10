/**
   @file contract_O1O2.c
   @brief perform the Diquarks - Baryon-Meson contraction
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_lr()
#include "gammas.h"             // gt_Gdag_gt()
#include "Ospinor.h"            // spinor_to_Ospinor()
#include "penta_contractions.h" // idx()
#include "spinmatrix_ops.h"     // spinmatrix_multiply()
#include "spinor_ops.h"         // spinmul_atomic_left()

// precompute the F-tensor, this is BY FAR the hottest part of this code
// any good ideas here would be more than welcome!
static void
precompute_F_O1O2_v2( double complex **F ,
		      const struct Ospinor OU1 ,
		      const struct Ospinor OD ,
		      const struct Ospinor OU2 ,
		      const struct Ospinor OM ,
		      const uint8_t **loc )
{
  // precompute all sub-spinmatrices
  struct spinmatrix temp5[ 9 ][ 9 ] , temp6[ 9 ][ 9 ] ;

  size_t c1 , c2 , c3 , c4 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      const struct spinmatrix t1 = OU1.C[ c1 ][ c2 ] ;
      const struct spinmatrix t3 = OU2.C[ c1 ][ c2 ] ; 
      for( c3 = 0 ; c3 < NC ; c3++ ) {
        for( c4 = 0 ; c4 < NC ; c4++ ) {
	  const struct spinmatrix t2 = OD.C[ c3 ][ c4 ] ;
	  const struct spinmatrix t4 = OM.C[ c3 ][ c4 ] ;
	  spinmatrix_multiply( temp5[c2+NC*c1][c4+NC*c3].D , t1.D , t2.D ) ;
	  spinmatrix_multiply( temp6[c2+NC*c1][c4+NC*c3].D , t3.D , t4.D ) ;
	  transpose_spinmatrix( temp6[c2+NC*c1][c4+NC*c3].D ) ;
	}
      }
    }
  }

  size_t i ;
  for( i = 0 ; i < PENTA_NCOLORS ; i++ ) {
        
    // preset the spinmatrices
    struct spinmatrix
      temp7 __attribute__ ((aligned(SPINT_ALIGNMENT))) ,
      temp8 __attribute__ ((aligned(SPINT_ALIGNMENT))) ;
    
    spinmatrix_multiply_T_avx( temp7.D ,
			       temp5[ loc[i][1] + NC*loc[i][0] ][ loc[i][3] + NC*loc[i][2] ].D ,
			       temp6[ loc[i][5] + NC*loc[i][4] ][ loc[i][7] + NC*loc[i][6] ].D ) ;
      
    spinmatrix_multiply_T_avx( temp8.D ,
			       temp5[ loc[i][5] + NC*loc[i][0] ][ loc[i][3] + NC*loc[i][2] ].D ,
			       temp6[ loc[i][1] + NC*loc[i][4] ][ loc[i][7] + NC*loc[i][6] ].D ) ;

    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + NS*d1 ][i] = temp7.D[d1][d2] - temp8.D[d1][d2] ;
      }
    }
  }

  return ;
}

// contract the colors of the F-tensor
static double complex
contract_colors_O1O2( const double complex *F )
{
  register double complex sum = 0.0 ;
  size_t b , c , g , h ; 
  for( b = 0 ; b < NC ; b++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      for( g = 0 ; g < NC ; g++ ) {
	for( h = 0 ; h < NC ; h++ ) {
	  // first set of epsilon identities
	  sum += ( +F[ idx( c , b , c , h , g , g , h , b ) ] 
		   -F[ idx( c , b , c , g , h , g , h , b ) ]
		   -F[ idx( g , b , c , h , c , g , h , b ) ]
		   +F[ idx( h , b , c , g , c , g , h , b ) ]
		   +F[ idx( g , b , c , c , h , g , h , b ) ]
		   -F[ idx( h , b , c , c , g , g , h , b ) ] ) ;
	  // second set of epsilon identities
	  sum -= ( +F[ idx( b , b , c , h , g , g , h , c ) ]
		   -F[ idx( b , b , c , g , h , g , h , c ) ]
		   -F[ idx( g , b , c , h , b , g , h , c ) ]
		   +F[ idx( h , b , c , g , b , g , h , c ) ]
		   +F[ idx( g , b , c , b , h , g , h , c ) ]
		   -F[ idx( h , b , c , b , g , g , h , c ) ] ) ;
	}
      }
    }
  }
  return sum ;
}

// 
void
contract_O1O2( struct spinmatrix *P ,
	       double complex **F ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS ,
	       const uint8_t **loc )
{
  // precompute some gammas
  struct gamma C1 = CGmu( OP1 , GAMMAS ) ;
  struct gamma C2 = CGmu( OP2 , GAMMAS ) ;
  struct gamma t2t = gt_Gdag_gt( OP2 , GAMMAS[ GAMMA_T ] ) ;
  struct gamma tC2t = gt_Gdag_gt( C2 , GAMMAS[ GAMMA_T ] ) ;
  
  // compute the common spinor "M" is like a b_s meson kinda
  struct spinor Temp = S ;
  struct spinor M = B ;
  gamma_mul_lr( &Temp , C1 , t2t ) ;
  spinmul_atomic_left( &M , Temp ) ;

  // precompute CG5 D \tilde{CG5}
  Temp = D ;
  gamma_mul_lr( &Temp , C1 , tC2t ) ;

  // convert to cache-friendly Ospinors
  struct Ospinor OU1 = spinor_to_Ospinor( U ) ;
  struct Ospinor OD = spinor_to_Ospinor( Temp ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( transpose_spinor( U ) ) ;
  struct Ospinor OM = spinor_to_Ospinor( M ) ;

  precompute_F_O1O2_v2( F , OU1 , OD , OU2 , OM , loc ) ;
  
  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {      
      P -> D[d1][d2] = contract_colors_O1O2( F[ d2 + NS*d1 ] ) ;
    }
  }
  
  return ;
}
