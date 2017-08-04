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
static double complex*
precompute_F_O1O1( const struct Ospinor OU1 ,
		   const struct Ospinor OU2 ,
		   const struct Ospinor OD ,
		   const struct Ospinor OS )
{
  double complex *F = malloc( 6561 * sizeof( double complex ) ) ;
  size_t i ;
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

    // this is the first contraction
    // T_S[ U_{bb'} G1 D_{cc'} ~G1 ] T_S[ U_{gg'} G2 S_{hh'} ~G2 ]
    struct spinmatrix
      temp1 = OU1.C[ idx[0] ][ idx[1] ] ,
      temp2 =  OD.C[ idx[2] ][ idx[3] ] ,
      temp3 = OU2.C[ idx[4] ][ idx[5] ] ,
      temp4 =  OS.C[ idx[6] ][ idx[7] ] ;
    
    F[i] =
      trace_prod_spinmatrices( temp1.D , temp2.D ) *
      trace_prod_spinmatrices( temp3.D , temp4.D ) ;

    // T_S[ U_{bg'} G1 D_{cc'} ~G1 U_{gb'} G2 S_{hh'} ~G2 ]
    temp1 = OU1.C[ idx[0] ][ idx[5] ] ;
    temp2 =  OD.C[ idx[2] ][ idx[3] ] ;
    temp3 = OU2.C[ idx[4] ][ idx[1] ] ;
    temp4 =  OS.C[ idx[6] ][ idx[7] ] ;

    // multiplty all of these together
    struct spinmatrix temp5 , temp6 ;
    spinmatrix_multiply( temp5.D , temp1.D , temp2.D ) ;
    spinmatrix_multiply( temp6.D , temp3.D , temp4.D ) ;
    
    F[i] -= trace_prod_spinmatrices( temp5.D , temp6.D ) ;
  }
  
  return F ;
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
  for( b = 0 ; b < NC ; b++ ) {
    for( prime = 0 ; prime < NC ; prime++ ) {   
      for( c = 0 ; c < NC ; c++ ) {
	for( g = 0 ; g < NC ; g++ ) {
	  for( h = 0 ; h < NC ; h++ ) {
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
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor bwdH ,
	       const struct gamma *GAMMAS )
{
  // Idea:: create a huge array with all the possible color components
  struct gamma CG5 = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  struct gamma tCG5t = gt_Gdag_gt( CG5 , GAMMAS[ GAMMA_T ] ) ;

  // temporary spinors
  struct spinor U1 = transpose_spinor( U ) ,
    U2 = transpose_spinor( U ) ,
    Dt = D , St = S , Bt = bwdH ;

  // perform some gamma multiplications
  gamma_mul_r( &U1 , CG5 ) ;
  gamma_mul_r( &Dt , tCG5t ) ;
  gamma_mul_r( &U2 , CG5 ) ;
  gamma_mul_r( &St , tCG5t ) ;
  gamma_mul_l( &Bt , GAMMAS[ GAMMA_T ] ) ;

  // switch to a color - dirac structure for much better
  // cache coherence
  struct Ospinor OU1 = spinor_to_Ospinor( U1 ) ;
  struct Ospinor OU2 = spinor_to_Ospinor( U2 ) ;
  struct Ospinor OD  = spinor_to_Ospinor( Dt ) ;
  struct Ospinor OS  = spinor_to_Ospinor( St ) ;
  
  // pre-compute all possible indices
  double complex *F = precompute_F_O1O1( OU1 , OU2 , OD , OS ) ;
  
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] = contract_colors( F , Bt.D[d1][d2] ) ;
    }
  } 
  free( F ) ;
  return ;
}
