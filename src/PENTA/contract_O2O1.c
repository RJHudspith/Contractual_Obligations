/**
   @file contract_O2O1.c
   @brief contract the Baryon-Meson with the Diquarks
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_lr()
#include "gammas.h"             // gt_Gdag_gt()
#include "Ospinor.h"            // spinor_to_Ospinor()
#include "penta_contractions.h" // idx()
#include "spinmatrix_ops.h"     // spinmatrix_multiply()
#include "spinor_ops.h"         // spinmul_atomic_left()

// precompute the F-tensor
static void
precompute_F_O2O1_v2( double complex **F ,
		      const struct Ospinor OM ,
		      const struct Ospinor OU1 ,
		      const struct Ospinor OD ,
		      const struct Ospinor OU2 ,
		      const uint8_t **loc )
{
  // precompute the spinmatrix sub-matrices
  struct spinmatrix temp5[ 9 ][ 9 ] , temp6[ 9 ][ 9 ] ;
  size_t c1 , c2 , c3 , c4 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      const struct spinmatrix t1 = OM.C[ c1 ][ c2 ] ;
      const struct spinmatrix t3 = OD.C[ c1 ][ c2 ] ;
      for( c3 = 0 ; c3 < NC ; c3++ ) {
        for( c4 = 0 ; c4 < NC ; c4++ ) {
	  const struct spinmatrix t2 = OU1.C[ c3 ][ c4 ] ;
	  const struct spinmatrix t4 = OU2.C[ c3 ][ c4 ] ;
	  spinmatrix_multiply( temp5[c2+NC*c1][c4+NC*c3].D , t1.D , t2.D ) ;  
	  spinmatrix_multiply( temp6[c2+NC*c1][c4+NC*c3].D , t3.D , t4.D ) ;
	  transpose_spinmatrix( temp6[c2+NC*c1][c4+NC*c3].D ) ;
	}
      }
    }
  }

  size_t i ;
  for( i = 0 ; i < PENTA_NCOLORS ; i++ ) {

    // there is a reproduction of results here that we can utilise
    // as this one is computed with swapped indices in a later evaluation
    struct spinmatrix \
      temp7 __attribute__((aligned(SPINT_ALIGNMENT))) ,
      temp8 __attribute__((aligned(SPINT_ALIGNMENT))) ;
    
    spinmatrix_multiply_T_avx( temp7.D ,
			       temp5[ loc[i][1] + NC*loc[i][0] ][ loc[i][3] + NC*loc[i][2] ].D ,
			       temp6[ loc[i][5] + NC*loc[i][4] ][ loc[i][7] + NC*loc[i][6] ].D ) ;

    spinmatrix_multiply_T_avx( temp8.D ,
			       temp5[ loc[i][1] + NC*loc[i][0] ][ loc[i][7] + NC*loc[i][2] ].D ,
			       temp6[ loc[i][5] + NC*loc[i][4] ][ loc[i][3] + NC*loc[i][6] ].D ) ;
    
    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] = temp7.D[d1][d2] - temp8.D[d1][d2] ;
      }
    }
  }

  return ;
}

// contract the colors of the F-tensor
static double complex
contract_colors_O2O1( const double complex *F )
{
  register double complex sum = 0.0 ;
  size_t a , b , c , prime ;
  for( a = 0 ; a < NC ; a++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      for( b = 0 ; b < NC ; b++ ) {
	for( prime = 0 ; prime < NC ; prime++ ) {
	  // first set of epsilon identities
	  sum += ( +F[ idx( prime , c , b , b , c , a , prime , a ) ]
		   -F[ idx( prime , b , c , b , c , a , prime , a ) ]
		   -F[ idx( prime , c , a , b , c , b , prime , a ) ]
		   +F[ idx( prime , b , a , b , c , c , prime , a ) ]
		   +F[ idx( prime , a , c , b , c , b , prime , a ) ]
		   -F[ idx( prime , a , b , b , c , c , prime , a ) ] ) ;
	  // second set of epsilon identities
	  sum -= ( +F[ idx( prime , c , b , b , c , prime , a , a ) ]
		   -F[ idx( prime , b , c , b , c , prime , a , a ) ]
		   -F[ idx( prime , c , a , b , c , prime , b , a ) ]
		   +F[ idx( prime , b , a , b , c , prime , c , a ) ]
		   +F[ idx( prime , a , c , b , c , prime , b , a ) ]
		   -F[ idx( prime , a , b , b , c , prime , c , a ) ] ) ;
	}
      }
    }
  }
  return sum ;
}

// 
void
contract_O2O1( struct spinmatrix *P ,
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
  struct gamma G5 = OP1 ;
  struct gamma C1 = CGmu( OP1 , GAMMAS ) ;

  struct gamma C2 = CGmu( OP2 , GAMMAS ) ;
  struct gamma tC2t = gt_Gdag_gt( C2 , GAMMAS[ GAMMA_T ] ) ;
  
  // compute the common spinor "M" is like the meson again
  struct spinor M = S ;
  gamma_mul_lr( &M , G5 , tC2t ) ;
  spinmul_atomic_left( &M , B ) ;
  //gamma_mul_l( &M , GAMMAS[ GAMMA_T ] ) ;
  
  // precompute
  struct spinor Temp = D ;
  gamma_mul_lr( &Temp , C1 , tC2t ) ;

  struct Ospinor OM  = spinor_to_Ospinor( M ) ;
  struct Ospinor OU1 = spinor_to_Ospinor( transpose_spinor( U ) ) ;
  struct Ospinor OD  = spinor_to_Ospinor( Temp ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( U ) ;

  // precompute the F-tensor
  precompute_F_O2O1_v2( F , OM , OU1 , OD , OU2 , loc ) ;
  
  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors_O2O1( F[ d2 + NS*d1 ] ) ;
    }
  }
  
  return ;
}
