/**
   @file bar_ops.c
   @brief baryon operations
 */
#include "common.h"

#include "spinor_ops.h" // spinor_zero_site()

#ifndef HAVE_EMMINTRIN_H

// This contracts the diquark with the remaining propagator
// This does the color trace Tr[ A . B ] 
const double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const int d0 ,
		 const int d1 ,
		 const int d2 ,
		 const int d3 )
{
  int c1, c2 ;
  register double corrr = 0.0 , corri = 0.0 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      corrr += creal( DiQ.D[d1][d0].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] ) 
	- cimag( DiQ.D[d1][d0].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] );
      corri += creal( DiQ.D[d1][d0].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] ) 
	+ cimag( DiQ.D[d1][d0].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] );
    }
  }
  return corrr + I * corri;
}

// cross color 
static inline void
cross_color( double complex *__restrict a ,
	     const double complex *__restrict b ,
	     const double complex *__restrict c ,
	     const int id1 ,
	     const int id2 )
{
#if NC == 3
  *a += b[ id1 + 0 ] * c[ id2 + 0 ] - b[ id2 + 0 ] * c[ id1 + 0 ] ; a++ ;
  *a += b[ id1 + 0 ] * c[ id2 + 3 ] - b[ id2 + 0 ] * c[ id1 + 3 ] ; a++ ;
  *a += b[ id1 + 0 ] * c[ id2 + 6 ] - b[ id2 + 0 ] * c[ id1 + 6 ] ; a++ ;

  *a += b[ id1 + 3 ] * c[ id2 + 0 ] - b[ id2 + 3 ] * c[ id1 + 0 ] ; a++ ;
  *a += b[ id1 + 3 ] * c[ id2 + 3 ] - b[ id2 + 3 ] * c[ id1 + 3 ] ; a++ ;
  *a += b[ id1 + 3 ] * c[ id2 + 6 ] - b[ id2 + 3 ] * c[ id1 + 6 ] ; a++ ;

  *a += b[ id1 + 6 ] * c[ id2 + 0 ] - b[ id2 + 6 ] * c[ id1 + 0 ] ; a++ ;
  *a += b[ id1 + 6 ] * c[ id2 + 3 ] - b[ id2 + 6 ] * c[ id1 + 3 ] ; a++ ;
  *a += b[ id1 + 6 ] * c[ id2 + 6 ] - b[ id2 + 6 ] * c[ id1 + 6 ] ; a++ ;
#else
  int c1 , c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      a[ c2 + c1 * NC ] += b[ id1 + c1 * NC ] * c[ id2 + NC * c2 ] ; 
      a[ c2 + c1 * NC ] -= b[ id2 + c1 * NC ] * c[ id1 + NC * c2 ] ;
    }
  }
#endif
}

// This carries out the color cross product and traces one set of Dirac indices.
// The result forms a diquark-type object
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S )
{
  // temporary 3-spinor space
  struct spinor T3SNK[ 3 ] ;

  // initialise temporary storage
  spinor_zero_site( &T3SNK[ 0 ] ) ;
  spinor_zero_site( &T3SNK[ 1 ] ) ;
  spinor_zero_site( &T3SNK[ 2 ] ) ;

  // Sink cross color and trace, this leaves only one set of dirac indices
  int i , j , d ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {

      double complex *a0 = (double complex*)T3SNK[0].D[i][j].C ;
      double complex *a1 = (double complex*)T3SNK[1].D[i][j].C ;
      double complex *a2 = (double complex*)T3SNK[2].D[i][j].C ;

      // Dirac sink trace
      for( d = 0 ; d < NS ; d++ ) {

	const double complex *b = (const double complex*)S.D[i][d].C ;
	const double complex *c = (const double complex*)DiQ->D[j][d].C ;

	cross_color( a0 , b , c , 1 , 2 ) ;
	cross_color( a1 , b , c , 2 , 0 ) ;
	cross_color( a2 , b , c , 0 , 1 ) ; 
      }
    }
  }

  // poke back into the diquark
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      DiQ -> D[i][j].C[0][0] = T3SNK[0].D[i][j].C[1][2] - T3SNK[0].D[i][j].C[2][1] ; 
      DiQ -> D[i][j].C[0][1] = T3SNK[1].D[i][j].C[1][2] - T3SNK[1].D[i][j].C[2][1] ; 
      DiQ -> D[i][j].C[0][2] = T3SNK[2].D[i][j].C[1][2] - T3SNK[2].D[i][j].C[2][1] ; 

      DiQ -> D[i][j].C[1][0] = T3SNK[0].D[i][j].C[2][0] - T3SNK[0].D[i][j].C[0][2] ; 
      DiQ -> D[i][j].C[1][1] = T3SNK[1].D[i][j].C[2][0] - T3SNK[1].D[i][j].C[0][2] ; 
      DiQ -> D[i][j].C[1][2] = T3SNK[2].D[i][j].C[2][0] - T3SNK[2].D[i][j].C[0][2] ; 

      DiQ -> D[i][j].C[2][0] = T3SNK[0].D[i][j].C[0][1] - T3SNK[0].D[i][j].C[1][0] ; 
      DiQ -> D[i][j].C[2][1] = T3SNK[1].D[i][j].C[0][1] - T3SNK[1].D[i][j].C[1][0] ;
      DiQ -> D[i][j].C[2][2] = T3SNK[2].D[i][j].C[0][1] - T3SNK[2].D[i][j].C[1][0] ;  
    }
  }

  return ;
}

#endif
