/**
   @file contract_O2O3.c
   @brief perform the contraction of our 2 baryon-meson operators
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
precompute_F_O2O3( const struct Ospinor OU1 ,
		   const struct Ospinor OM ,
		   const struct Ospinor OU2 )
{
  double complex **F = malloc( NSNS * sizeof( double complex* ) ) ;
  size_t i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    F[i] = malloc( 729 * sizeof( double complex ) ) ;
  }

  for( i = 0 ; i < 729 ; i++ ) {
    
    // map the indices
    size_t j , idx[ 6 ] , sub = NC , div = 1 ;
    for( j = 6 ; j > 0 ; j-- ) {
      idx[ 6 - j ] = ( i % sub ) / div ;
      sub *= NC ;
      div *= NC ;
    }
    
    // preset the spinmatrices
    struct spinmatrix
      temp1 = OU1.C[ idx[0] ][ idx[1] ] ,
      temp2 = OM.C[ idx[2] ][ idx[3] ] ,
      temp3 = OU2.C[ idx[4] ][ idx[5] ] ,
      temp4 , temp5 ;

    const double complex tr = trace_prod_spinmatrices( temp1.D , temp2.D ) ;

    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] = tr * temp3.D[d1][d2] ;
      }
    }

    spinmatrix_multiply( temp4.D , temp1.D , temp2.D ) ;
    spinmatrix_multiply( temp5.D , temp4.D , temp3.D ) ;

    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] = temp5.D[d2][d1] ;
      }
    }
  }
  return F ;
}


// contract the colors of the F-tensor
static double complex
contract_colors_O2O3( const double complex *F )
{
  register double complex sum = 0.0 ;
  size_t a , b , c ;
  for( a = 0 ; a < NC ; a++ ) {
    for( b = 0 ; b < NC ; b++ ) {
      for( c = 0 ; c < NC ; c++ ) {
	// first set of epsilon identities
	sum += ( +F[ idx2( b , b , c , c , a , a ) ] 
		 -F[ idx2( c , b , c , b , a , a ) ]
		 -F[ idx2( a , b , c , c , a , b ) ]
		 +F[ idx2( a , b , c , b , a , c ) ]
		 +F[ idx2( c , b , c , a , a , b ) ]
		 -F[ idx2( b , b , c , a , a , c ) ] ) ;
      }
    }
  }

  return sum ;
}
// 
void
contract_O2O3( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  // precompute some gammas
  struct gamma G5 = GAMMAS[ GAMMA_5 ] ;
  struct gamma tG5t = gt_Gdag_gt( G5 , GAMMAS[ GAMMA_T ] ) ;
  
  struct gamma CG5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma tCG5t = gt_Gdag_gt( CG5 , GAMMAS[ GAMMA_T ] ) ;

  // precompute [ (Cg5 D tg5t) B (g5 S tCg5t) ]
  struct spinor M = S ;
  gamma_mul_lr( &M , G5 , tCG5t ) ;
  struct spinor temp = D ;
  gamma_mul_lr( &temp , CG5 , tG5t ) ;
  
  spinmul_atomic_left( &M , B ) ;
  spinmul_atomic_left( &M , temp ) ;

  const struct Ospinor OU1 = spinor_to_Ospinor( U ) ;
  const struct Ospinor OU2 = spinor_to_Ospinor( U ) ;
  const struct Ospinor OM = spinor_to_Ospinor( U ) ;
  
  double complex **F = precompute_F_O2O3( OU1 , OM , OU2 ) ;

  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors_O2O3( F[ d2 + NS*d1 ] ) ;
    }
  }
  
  // free the F-tensor
  for( d1 = 0 ; d1 < NSNS ; d1++ ) {
    free( F[ d1 ] ) ;
  }
  free( F ) ;
  return ;
}
