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

// precompute the F-tensor
static double complex**
precompute_F_O1O2( const struct Ospinor OU1 ,
		   const struct Ospinor OD ,
		   const struct Ospinor OU2 ,
		   const struct Ospinor OM )
{
  double complex **F = malloc( NSNS * sizeof( double complex* ) ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    F[i] = malloc( 6561 * sizeof( double complex ) ) ;
  }
    
  for( i = 0 ; i < 6561 ; i++ ) {

    // map the indices correctly
    // b = idx[0] , b' = idx[1]
    // c = idx[2] , c' = idx[3]
    // g = idx[4] , g' = idx[5]
    // h = idx[6] , h' = idx[7]
    size_t j , idx[ 8 ] , sub = NC , div = 1 ;
    for( j = 8 ; j > 0 ; j-- ) {
      idx[ 8 - j ] = ( i % sub ) / div ;
      sub *= NC ;
      div *= NC ;
    }

    // preset the spinmatrices
    struct spinmatrix
      temp1 = OU1.C[ idx[0] ][ idx[1] ] ,
      temp2 = OD.C[ idx[2] ][ idx[3] ] ,
      temp3 = OU2.C[ idx[4] ][ idx[5] ] ,
      temp4 = OM.C[ idx[6] ][ idx[7] ] ,
      temp5 , temp6 , temp7 ;

    spinmatrix_multiply( temp5.D , temp1.D , temp2.D ) ;
    spinmatrix_multiply( temp6.D , temp3.D , temp4.D ) ;
    spinmatrix_multiply( temp7.D , temp5.D , temp6.D ) ;
    
    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] = temp7.D[d1][d2] ;
      }
    }

    temp1 = OU1.C[ idx[0] ][ idx[5] ] ;
    temp3 = OU2.C[ idx[4] ][ idx[1] ] ;

    spinmatrix_multiply( temp5.D , temp1.D , temp2.D ) ;
    spinmatrix_multiply( temp6.D , temp3.D , temp4.D ) ;
    spinmatrix_multiply( temp7.D , temp5.D , temp6.D ) ;

    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] -= temp7.D[d1][d2] ;
      }
    }
  }
  
  return F ;
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
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  // precompute some gammas
  struct gamma CG5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma tG5t = gt_Gdag_gt( GAMMAS[ GAMMA_5 ] , GAMMAS[ GAMMA_T ] ) ;
  struct gamma tCG5t = gt_Gdag_gt( CG5 , GAMMAS[ GAMMA_T ] ) ;
  
  // compute the common spinor "M" is like a b_s meson kinda
  struct spinor Temp = S ;
  struct spinor M = B ;
  gamma_mul_lr( &Temp , CG5 , tG5t ) ;
  spinmul_atomic_left( &M , Temp ) ;

  // precompute CG5 D \tilde{CG5}
  Temp = D ;
  gamma_mul_lr( &Temp , CG5 , tCG5t ) ;

  // convert to cache-friendly Ospinors
  struct Ospinor OU1 = spinor_to_Ospinor( U ) ;
  struct Ospinor OD = spinor_to_Ospinor( Temp ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( transpose_spinor( U ) ) ;
  struct Ospinor OM = spinor_to_Ospinor( M ) ;

  double complex **F = precompute_F_O1O2( OU1 , OD , OU2 , OM ) ;

  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors_O1O2( F[ d2 + NS*d1 ] ) ;
    }
  }

  // free the F-tensor
  for( d1 = 0 ; d1 < NSNS ; d1++ ) {
    free( F[ d1 ] ) ;
  }
  free( F ) ;
  
  return ;
}
