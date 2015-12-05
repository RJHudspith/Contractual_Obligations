/**
   @file basis_conversions.c
   @brief conversions between gamma bases

   TODO :: cache coherency - J
 */

#include "common.h"

// convert chiral to Nrel -- updated ordering
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  size_t c1c2 ;
  for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
    const size_t c1 = c1c2 / NC ;
    const size_t c2 = c1c2 % NC ;
    // ok, I get something different than we had before
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

// work in progress - J
#if 0 

typedef double complex dcl ;

// okey dokey there are sixteen possible combinations of plusses and minuses
// encoded as binary, plus is 0 minus is 1
const dcl  f0( dcl a , dcl b , dcl c , dcl d ) { return  a+b+c+d ; } // ++++
const dcl  f1( dcl a , dcl b , dcl c , dcl d ) { return  a+b+c-d ; } // +++-
const dcl  f2( dcl a , dcl b , dcl c , dcl d ) { return  a+b-c+d ; } // ++-+
const dcl  f3( dcl a , dcl b , dcl c , dcl d ) { return  a+b-c-d ; } // ++--
const dcl  f4( dcl a , dcl b , dcl c , dcl d ) { return  a-b+c+d ; } // +-++
const dcl  f5( dcl a , dcl b , dcl c , dcl d ) { return  a-b+c-d ; } // +-+-
const dcl  f6( dcl a , dcl b , dcl c , dcl d ) { return  a-b-c+d ; } // +--+
const dcl  f7( dcl a , dcl b , dcl c , dcl d ) { return  a-b-c-d ; } // +---
const dcl  f8( dcl a , dcl b , dcl c , dcl d ) { return -a+b+c+d ; } // -+++
const dcl  f9( dcl a , dcl b , dcl c , dcl d ) { return -a+b+c-d ; } // -++-
const dcl f10( dcl a , dcl b , dcl c , dcl d ) { return -a+b-c+d ; } // -+-+
const dcl f11( dcl a , dcl b , dcl c , dcl d ) { return -a+b-c-d ; } // -+--
const dcl f12( dcl a , dcl b , dcl c , dcl d ) { return -a-b+c+d ; } // --++
const dcl f13( dcl a , dcl b , dcl c , dcl d ) { return -a-b+c-d ; } // --+-
const dcl f14( dcl a , dcl b , dcl c , dcl d ) { return -a-b-c+d ; } // ---+
const dcl f15( dcl a , dcl b , dcl c , dcl d ) { return -a-b-c+d ; } // ----

// actual matrix
static void
rotate_dmatrix( dcl *res ,
		const dcl *P1 , 
		const dcl *P2 ,
		const dcl *P3 ,
		const dcl *P4 ,
		const dcl (*f)( dcl a , 
				dcl b , 
				dcl c , 
				dcl d ) )
{
  size_t c1c2 ;
  for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
    *res = 0.5 * f( *P1 , *P2 , *P3 , *P4 ) ;
    res++ , P1++ , P2++ , P3++ , P4++ ;
  }
  return ;
}
 

// convert chiral to non-relativistic
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  // first row
  rotate_dmatrix( (dcl*)S -> D[0][0].C , (dcl*)P.D[1][1].C , (dcl*)P.D[1][3].C , (dcl*)P.D[3][1].C , (dcl*)P.D[3][3].C , f6  ) ;
  rotate_dmatrix( (dcl*)S -> D[0][1].C , (dcl*)P.D[1][0].C , (dcl*)P.D[1][2].C , (dcl*)P.D[3][0].C , (dcl*)P.D[3][2].C , f6  ) ;
  rotate_dmatrix( (dcl*)S -> D[0][2].C , (dcl*)P.D[1][1].C , (dcl*)P.D[1][3].C , (dcl*)P.D[3][1].C , (dcl*)P.D[3][3].C , f12 ) ;
  rotate_dmatrix( (dcl*)S -> D[0][3].C , (dcl*)P.D[1][0].C , (dcl*)P.D[1][2].C , (dcl*)P.D[3][0].C , (dcl*)P.D[3][2].C , f12 ) ;
  // second row
  rotate_dmatrix( (dcl*)S -> D[1][0].C , (dcl*)P.D[0][1].C , (dcl*)P.D[0][3].C , (dcl*)P.D[2][1].C , (dcl*)P.D[2][3].C , f6  ) ;
  rotate_dmatrix( (dcl*)S -> D[1][1].C , (dcl*)P.D[0][0].C , (dcl*)P.D[0][2].C , (dcl*)P.D[2][0].C , (dcl*)P.D[2][2].C , f6  ) ;
  rotate_dmatrix( (dcl*)S -> D[1][2].C , (dcl*)P.D[0][1].C , (dcl*)P.D[0][3].C , (dcl*)P.D[2][1].C , (dcl*)P.D[2][3].C , f12 ) ;
  rotate_dmatrix( (dcl*)S -> D[1][3].C , (dcl*)P.D[0][0].C , (dcl*)P.D[0][2].C , (dcl*)P.D[2][0].C , (dcl*)P.D[2][2].C , f12 ) ;
  // third row
  rotate_dmatrix( (dcl*)S -> D[2][0].C , (dcl*)P.D[1][1].C , (dcl*)P.D[1][3].C , (dcl*)P.D[3][1].C , (dcl*)P.D[3][3].C , f10 ) ;
  rotate_dmatrix( (dcl*)S -> D[2][1].C , (dcl*)P.D[1][0].C , (dcl*)P.D[1][2].C , (dcl*)P.D[3][0].C , (dcl*)P.D[3][2].C , f10 ) ;
  rotate_dmatrix( (dcl*)S -> D[2][2].C , (dcl*)P.D[1][1].C , (dcl*)P.D[1][3].C , (dcl*)P.D[3][1].C , (dcl*)P.D[3][3].C , f0  ) ;
  rotate_dmatrix( (dcl*)S -> D[2][3].C , (dcl*)P.D[1][0].C , (dcl*)P.D[1][2].C , (dcl*)P.D[3][0].C , (dcl*)P.D[3][2].C , f0  ) ;
  // fourth row
  rotate_dmatrix( (dcl*)S -> D[3][0].C , (dcl*)P.D[0][1].C , (dcl*)P.D[0][3].C , (dcl*)P.D[2][1].C , (dcl*)P.D[2][3].C , f10 ) ;
  rotate_dmatrix( (dcl*)S -> D[3][1].C , (dcl*)P.D[0][0].C , (dcl*)P.D[0][2].C , (dcl*)P.D[2][0].C , (dcl*)P.D[2][2].C , f10 ) ;
  rotate_dmatrix( (dcl*)S -> D[3][2].C , (dcl*)P.D[0][1].C , (dcl*)P.D[0][3].C , (dcl*)P.D[2][1].C , (dcl*)P.D[2][3].C , f0  ) ;
  rotate_dmatrix( (dcl*)S -> D[3][3].C , (dcl*)P.D[0][0].C , (dcl*)P.D[0][2].C , (dcl*)P.D[2][0].C , (dcl*)P.D[2][2].C , f0  ) ;
  return ;
}
#endif
