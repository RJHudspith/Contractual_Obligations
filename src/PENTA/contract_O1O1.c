/**
   @file contract_O1O1.c
   @brief perform contractions for the diquarks
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_r()
#include "gammas.h"             // Cgmu()
#include "Ospinor.h"            // opposite ordering spinor
#include "penta_contractions.h" // idx()
#include "spinmatrix_ops.h"     // trace_prod, get_spinmatrix()
#include "spinor_ops.h"         // transpose_spinor()

// precompute the F-tensor
static void
precompute_F_O1O1_v2( double complex *F ,
		      const struct Ospinor OU1 ,
		      const struct Ospinor OU2 ,
		      const struct Ospinor OD ,
		      const struct Ospinor OS ,
		      const uint8_t **loc )
{
  // precompute all sub-spinmatrices
  struct spinmatrix temp5[9][9] , temp6[9][9] ;
	    
  size_t c1 , c2 , c3 , c4 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      
      const struct spinmatrix t1 = OU1.C[ c1 ][ c2 ] ;
      const struct spinmatrix t3 = OU2.C[ c1 ][ c2 ] ;
      
      for( c3 = 0 ; c3 < NC ; c3++ ) {
        for( c4 = 0 ; c4 < NC ; c4++ ) {

	  const struct spinmatrix t2 = OD.C[ c3 ][ c4 ] ;
	  const struct spinmatrix t4 = OS.C[ c3 ][ c4 ] ;
	  
	  spinmatrix_multiply( temp5[c2+NC*c1][c4+NC*c3].D , t1.D , t2.D ) ;
	  
	  spinmatrix_multiply( temp6[c2+NC*c1][c4+NC*c3].D , t3.D , t4.D ) ;
	}
      }
    }
  }
  
  size_t i ;
  for( i = 0 ; i < PENTA_NCOLORS ; i++ ) {
    
    F[i] =
      spinmatrix_trace( temp5[ loc[i][1] + NC*loc[i][0] ][ loc[i][3] + NC*loc[i][2] ].D ) *
      spinmatrix_trace( temp6[ loc[i][5] + NC*loc[i][4] ][ loc[i][7] + NC*loc[i][6] ].D ) ;

    F[i] -=
      trace_prod_spinmatrices( temp5[ loc[i][5] + NC*loc[i][0] ][ loc[i][3] + NC*loc[i][2] ].D ,
			       temp6[ loc[i][1] + NC*loc[i][4] ][ loc[i][7] + NC*loc[i][6] ].D ) ;
  }
  
  return ;
}

// perform a color contraction
static double complex
contract_colors( const double complex *F ,
		 const struct colormatrix B )
{
  // it doesn't matter which order of indices you
  // sum over and there is always a primed one hitting the b
  // so I stick that out the front, logically there would be
  // 5 color sums for our pentaquark
  register double complex sum = 0.0 ;
  size_t b , c , g , h , prime ;
  for( h = 0 ; h < NC ; h++ ) {
    for( g = 0 ; g < NC ; g++ ) {
      for( c = 0 ; c < NC ; c++ ) {
	for( prime = 0 ; prime < NC ; prime++ ) {
	  for( b = 0 ; b < NC ; b++ ) {
	    // op1
	    sum +=
	      ( F[ idx( b , prime , c , c , g , g , h , h ) ] -
		F[ idx( b , prime , c , c , g , h , h , g ) ] -
		F[ idx( b , prime , c , g , g , c , h , h ) ] +
		F[ idx( b , prime , c , h , g , c , h , g ) ] +
		F[ idx( b , prime , c , g , g , h , h , c ) ] -
		F[ idx( b , prime , c , h , g , g , h , c ) ] ) *
	      B.C[ prime ][ b ] ;
	    // op2
	    sum -=
	      ( F[ idx( b , c , c , prime , g , g , h , h ) ] -
		F[ idx( b , c , c , prime , g , h , h , g ) ] -
		F[ idx( b , g , c , prime , g , c , h , h ) ] +
		F[ idx( b , h , c , prime , g , c , h , g ) ] +
		F[ idx( b , g , c , prime , g , h , h , c ) ] -
		F[ idx( b , h , c , prime , g , g , h , c ) ] ) *
	      B.C[ prime ][ b ] ;
	    // op3
	    sum -=
	      ( F[ idx( b , prime , c , b , g , g , h , h ) ] -
		F[ idx( b , prime , c , b , g , h , h , g ) ] -
		F[ idx( b , prime , c , g , g , b , h , h ) ] +
		F[ idx( b , prime , c , h , g , b , h , g ) ] +
		F[ idx( b , prime , c , g , g , h , h , b ) ] -
		F[ idx( b , prime , c , h , g , g , h , b ) ] ) *
	      B.C[ prime ][ c ] ;    
	    // op4
	    sum +=
	      ( F[ idx( b , b , c , prime , g , g , h , h ) ] -
		F[ idx( b , b , c , prime , g , h , h , g ) ] -
		F[ idx( b , g , c , prime , g , b , h , h ) ] +
		F[ idx( b , h , c , prime , g , b , h , g ) ] +
		F[ idx( b , g , c , prime , g , h , h , b ) ] -
		F[ idx( b , h , c , prime , g , g , h , b ) ] ) *
	      B.C[ prime ][ c ] ;

	  }
	}
      }
    }
  }
  return sum ;
}

// (ud)(us)\bar{b} with ud us diquarks
void
contract_O1O1( struct spinmatrix *P ,
	       double complex **F ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor bwdH ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS ,
	       const uint8_t **loc )
{
  // Idea:: create a huge array with all the possible color components
  struct gamma C1   = CGmu( OP1 , GAMMAS ) ;
  struct gamma tC1t = gt_Gdag_gt( C1 , GAMMAS[ GAMMA_T ] ) ;
  struct gamma C2   = CGmu( OP2 , GAMMAS ) ;
  struct gamma tC2t = gt_Gdag_gt( C2 , GAMMAS[ GAMMA_T ] ) ;

  // temporary spinors
  struct spinor U1 = transpose_spinor( U ) ,
    U2 = transpose_spinor( U ) ,
    Dt = D , St = S , Bt = bwdH ;

  // perform some gamma multiplications
  gamma_mul_r( &U1 , C1 ) ;
  gamma_mul_r( &Dt , tC1t ) ;
  gamma_mul_r( &U2 , C2 ) ;
  gamma_mul_r( &St , tC2t ) ;
  //gamma_mul_l( &Bt , GAMMAS[ GAMMA_T ] ) ;

  // switch to a color - dirac structure for much better
  // cache coherence
  struct Ospinor OU1 = spinor_to_Ospinor( U1 ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( U2 ) ;
  struct Ospinor OD  = spinor_to_Ospinor( Dt ) ;
  struct Ospinor OS  = spinor_to_Ospinor( St ) ;
  
  // pre-compute all possible indices
  precompute_F_O1O1_v2( F[0] , OU1 , OU2 , OD , OS , loc ) ;

  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors( F[0] , Bt.D[d1][d2] ) ;
    }
  } 

  return ;
}
