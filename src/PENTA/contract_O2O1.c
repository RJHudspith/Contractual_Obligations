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
static double complex**
precompute_F_O2O1( const struct Ospinor OM ,
		   const struct Ospinor OU1 ,
		   const struct Ospinor OD ,
		   const struct Ospinor OU2 )
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
      temp1 = OM.C[ idx[0] ][ idx[1] ] ,
      temp2 = OU1.C[ idx[2] ][ idx[3] ] ,
      temp3 = OD.C[ idx[4] ][ idx[5] ] ,
      temp4 = OU2.C[ idx[6] ][ idx[7] ] ,
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

    temp2 = OU1.C[ idx[2] ][ idx[7] ] ;
    temp4 = OU2.C[ idx[6] ][ idx[3] ] ;

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
contract_colors_O2O1( const double complex *F )
{
  register double complex sum = 0.0 ;
  size_t a , b , c , prime ;
  for( prime = 0 ; prime < NC ; prime++ ) {
    for( b = 0 ; b < NC ; b++ ) {
      for( c = 0 ; c < NC ; c++ ) {
	for( a = 0 ; a < NC ; a++ ) {
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
	       const struct spinor L ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  struct gamma G5 = GAMMAS[ GAMMA_5 ] ;
  struct gamma CG5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma tCG5t = gt_Gdag_gt( CG5 , GAMMAS[ GAMMA_T ] ) ;
  
  // compute the common spinor "M"
  struct spinor M = S ;
  gamma_mul_lr( &M , G5 , tCG5t ) ;
  spinmul_atomic_left( &M , B ) ;
  gamma_mul_l( &M , GAMMAS[ GAMMA_T ] ) ;
  
  // precompute
  struct spinor Temp = L ;
  gamma_mul_lr( &Temp , CG5 , tCG5t ) ;

  struct Ospinor OM  = spinor_to_Ospinor( M ) ;
  struct Ospinor OU1 = spinor_to_Ospinor( transpose_spinor( L ) ) ;
  struct Ospinor OD  = spinor_to_Ospinor( Temp ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( L ) ;

  double complex **F = precompute_F_O2O1( OM , OU1 , OD , OU2 ) ;

  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors_O2O1( F[ d2 + NS*d1 ] ) ;
    }
  }

  // free the F-tensor
  for( d1 = 0 ; d1 < NSNS ; d1++ ) {
    free( F[ d1 ] ) ;
  }
  free( F ) ;
  
  return ;
}
