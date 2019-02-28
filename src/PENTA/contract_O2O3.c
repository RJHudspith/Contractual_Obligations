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
static void
precompute_F_O2O3_v2( double complex **F ,
		      const struct Ospinor OUT ,
		      const struct Ospinor OM ,
		      const struct Ospinor OU ,
		      const uint8_t **loc )
{
  struct spinmatrix temp4[ 9 ][ 9 ] ;
  size_t c1 , c2 , c3 , c4 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      struct spinmatrix t1 = OUT.C[ c1 ][ c2 ] ;
      for( c3 = 0 ; c3 < NC ; c3++ ) {
	for( c4 = 0 ; c4 < NC ; c4++ ) {
	  struct spinmatrix t2 = OM.C[ c3 ][ c4 ] ;
	  spinmatrix_multiply( temp4[c2+c1*NC][c4+c3*NC].D ,
			       t1.D , t2.D ) ;
	}
      }
    }
  }

  // loop possible index combinations
  size_t i ;
  for( i = 0 ; i < 729 ; i++ ) {
      
    // preset the spinmatrices
    struct spinmatrix temp3 = OU.C[ loc[i][4] ][ loc[i][5] ] ;
    struct spinmatrix temp5 ;

    size_t c1 = loc[i][0]+loc[i][1]*NC ;
    size_t c2 = loc[i][2]+loc[i][3]*NC ;
    const double complex tr = spinmatrix_trace( temp4[c1][c2].D ) ;

    size_t d1 , d2 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] = tr * temp3.D[d1][d2] ;
      }
    }

    // this does the cross term NOTE :: both U props are transposed
    temp3 = OUT.C[ loc[i][0] ][ loc[i][4] ] ;
    c1 = loc[i][5]+loc[i][1]*NC ;
    
    spinmatrix_multiply( temp5.D , temp4[c1][c2].D , temp3.D ) ;

    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	F[ d2 + d1*NS ][i] -= temp5.D[d2][d1] ;
      }
    }
  }

  return ;
}

// contract the colors of the F-tensor
static double complex
contract_colors_O2O3( const double complex *F )
{
  register double complex sum = 0.0 ;
  size_t a , b , c ;
  for( a = 0 ; a < NC ; a++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      for( b = 0 ; b < NC ; b++ ) {
	// epsilon identities
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
  // here OP1 is the gamma for the both the diquark and meson
  // and OP2 is the gamma structure for the daggered guy
  struct gamma G1    = OP1 ;
  struct gamma CG1   = CGmu( OP1 , GAMMAS ) ;

  struct gamma G2    = OP2 ;
  struct gamma tG2t  = gt_Gdag_gt( G2 , GAMMAS[ GAMMA_T ] ) ;
  struct gamma tCG2t = gt_Gdag_gt( CGmu( G2 , GAMMAS ) ,
				   GAMMAS[ GAMMA_T ] ) ;
  
  // precompute [ (CG1 D tG2t) B (G1 S tCG2t) ]

  // LHS of the B
  struct spinor temp = D ;
  gamma_mul_lr( &temp , CG1 , tG2t ) ;

  // RHS of the B
  struct spinor M = S ;
  gamma_mul_lr( &M , G1 , tCG2t ) ;
 
  spinmul_atomic_left( &M , B ) ;
  spinmul_atomic_left( &M , temp ) ;

  const struct Ospinor OUT = spinor_to_Ospinor( transpose_spinor( U ) ) ;
  const struct Ospinor OU  = spinor_to_Ospinor( U ) ;
  const struct Ospinor OM  = spinor_to_Ospinor( M ) ;

  precompute_F_O2O3_v2( F , OUT , OM , OU , loc ) ;

  // compute the color contraction
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors_O2O3( F[ d2 + NS*d1 ] ) ;
    }
  }
  
  return ;
}
