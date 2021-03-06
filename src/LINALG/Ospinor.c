/**
   @file Ospinor.c
   @brief Opposite-ordered spinor operations
 */
#include "common.h"

#include "spinmatrix_ops.h"

// copy a spinor to an opposite-indexed spinor
struct Ospinor
spinor_to_Ospinor( const struct spinor S )
{
  struct Ospinor O ;
  size_t d1 , d2 , c1 ,c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      for( d1 = 0 ; d1 < NS ; d1++ ) {
	for( d2 = 0 ; d2 < NS ; d2++ ) {
	  O.C[c1][c2].D[d1][d2] = S.D[d1][d2].C[c1][c2] ;
	}
      }
    }
  }
  return O ;
}

void
gamma_mul_r_Ospinor( struct Ospinor *S ,
		     const struct gamma Gamma )
{
  size_t d1 , d2 , c1 ,c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      spinmatrix_gamma( (void*)S -> C[c1][c2].D , Gamma ) ; 
    }
  }
  return ;
}

// multiplies by a gamma on the right and takes the traspose
void
gamma_mul_r_OspinorT( struct Ospinor *S ,
		      const struct gamma Gamma )
{
  size_t d1 , d2 , c1 ,c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      spinmatrix_gamma( (void*)S -> C[c1][c2].D , Gamma ) ;
      transpose_spinmatrix( (void*)S -> C[c1][c2].D ) ;
    }
  }
  return ;
}
