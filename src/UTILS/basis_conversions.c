/**
   @file basis_conversions.c
   @brief conversions between gamma bases

   TODO :: cache coherency - J
 */
#include "common.h"

// convert chiral to Nrel
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  size_t c1c2 ;
  for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
    const size_t c1 = c1c2 / NC ;
    const size_t c2 = c1c2 % NC ;
    /* G&L rotation
    // first row
    S -> D[0][0].C[c1][c2] = 0.5 * ( +P.D[1][1].C[c1][c2] - P.D[1][3].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[0][1].C[c1][c2] = 0.5 * ( +P.D[1][0].C[c1][c2] - P.D[1][2].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    S -> D[0][2].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[0][3].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    // second row
    S -> D[1][0].C[c1][c2] = 0.5 * ( +P.D[0][1].C[c1][c2] - P.D[0][3].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[1][1].C[c1][c2] = 0.5 * ( +P.D[0][0].C[c1][c2] - P.D[0][2].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    S -> D[1][2].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[1][3].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    // third row
    S -> D[2][0].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] + P.D[1][3].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[2][1].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] + P.D[1][2].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    S -> D[2][2].C[c1][c2] = 0.5 * ( +P.D[1][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[2][3].C[c1][c2] = 0.5 * ( +P.D[1][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    // final row
    S -> D[3][0].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] + P.D[0][3].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[3][1].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] + P.D[0][2].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    S -> D[3][2].C[c1][c2] = 0.5 * ( +P.D[0][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[3][3].C[c1][c2] = 0.5 * ( +P.D[0][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    */

   // first row
    S -> D[0][0].C[c1][c2] = 0.5 * ( +P.D[1][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[0][1].C[c1][c2] = 0.5 * ( +P.D[1][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    S -> D[0][2].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] + P.D[1][3].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[0][3].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] + P.D[1][2].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    // second row
    S -> D[1][0].C[c1][c2] = 0.5 * ( +P.D[0][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[1][1].C[c1][c2] = 0.5 * ( +P.D[0][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    S -> D[1][2].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] + P.D[0][3].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[1][3].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] + P.D[0][2].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    // third row
    S -> D[2][0].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[2][1].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    S -> D[2][2].C[c1][c2] = 0.5 * ( +P.D[1][1].C[c1][c2] - P.D[1][3].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[3][3].C[c1][c2] ) ;
    S -> D[2][3].C[c1][c2] = 0.5 * ( +P.D[1][0].C[c1][c2] - P.D[1][2].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[3][2].C[c1][c2] ) ;
    // final row
    S -> D[3][0].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[3][1].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;
    S -> D[3][2].C[c1][c2] = 0.5 * ( +P.D[0][1].C[c1][c2] - P.D[0][3].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[2][3].C[c1][c2] ) ;
    S -> D[3][3].C[c1][c2] = 0.5 * ( +P.D[0][0].C[c1][c2] - P.D[0][2].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[2][2].C[c1][c2] ) ;

  }
  return ;
}  

// rotate a timeslice
void
nrel_rotate_slice( struct spinor *S )
{
  size_t site ;
#pragma omp parallel for private(site) 
  for( site = 0 ; site < LCU ; site++ ) {
    chiral_to_nrel( &S[ site ] ) ;
  }
  return ;
}

// rotate depending on proptype
void
rotate_offdiag_2( struct spinor *S1 ,
		  const proptype basis1 ,
		  struct spinor *S2 ,
		  const proptype basis2 )
{
  // if we are doing nonrel-chiral mesons we switch chiral to nrel
  if( basis1 == CHIRAL && ( basis2 == NREL_FWD || 
			    basis2 == NREL_BWD ) ) {
    nrel_rotate_slice( S1 ) ;
  } else if( basis2 == CHIRAL && ( basis1 == NREL_FWD || 
				   basis1 == NREL_BWD ) ) {
    nrel_rotate_slice( S2 ) ;
  }
  return ;
}

// rotate depending on proptype
void
rotate_offdiag_3( struct spinor *S1 ,
		  const proptype basis1 ,
		  struct spinor *S2 ,
		  const proptype basis2 , 
		  struct spinor *S3 ,
		  const proptype basis3 )
{
  if( basis1 == CHIRAL && ( basis2 == NREL_FWD || basis2 == NREL_BWD || 
			    basis3 == NREL_FWD || basis3 == NREL_BWD ) ) {
    nrel_rotate_slice( S1 ) ;
  } 
  if( basis2 == CHIRAL && ( basis1 == NREL_FWD || basis1 == NREL_BWD || 
			    basis3 == NREL_FWD || basis3 == NREL_BWD ) ) {
    nrel_rotate_slice( S2 ) ;
  } 
  if( basis3 == CHIRAL && ( basis1 == NREL_FWD || basis1 == NREL_BWD ||
			    basis2 == NREL_FWD || basis2 == NREL_BWD ) ) {
    nrel_rotate_slice( S3 ) ;
  }
  return ;
}
