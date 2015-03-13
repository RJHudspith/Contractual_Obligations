/**
   @file basis_conversions.c
   @brief conversions between gamma bases
 */

#include "common.h"

// convert chiral to non-relativistic
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  int c1 , c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      S -> D[0][0].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[1][0].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
      S -> D[2][0].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[3][0].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

      S -> D[0][1].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[1][1].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
      S -> D[2][1].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[3][1].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );

      S -> D[0][2].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[1][2].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
      S -> D[2][2].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[3][2].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

      S -> D[0][3].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[1][3].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
      S -> D[2][3].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[3][3].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );	
    }
  }
  return ;
}

// rotate a timeslice
void
nrel_rotate_slice( struct spinor *S )
{
  int site ;
#pragma omp parallel for private(site) 
  for( site = 0 ; site < LCU ; site++ ) {
    chiral_to_nrel( &S[ site ] ) ;
  }
  return ;
}
