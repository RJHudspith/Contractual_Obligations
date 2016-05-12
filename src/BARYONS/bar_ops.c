/**
   @file bar_ops.c
   @brief baryon operations

   TODO :: try and generalise for SU(N)
 */
#include "common.h"

#include "spinor_ops.h" // spinor_zero_site()

#ifndef HAVE_EMMINTRIN_H

// This contracts the diquark with the remaining propagator
// This effectively does the color trace Tr[ A . B^T ] 
double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const size_t d0 ,
		 const size_t d1 ,
		 const size_t d2 ,
		 const size_t d3 )
{
  register double corrr = 0.0 , corri = 0.0 ;
  const double complex *diq = (const double complex*)DiQ.D[d0][d1].C ;
  const double complex *s = (const double complex*)S.D[d2][d3].C ;
  size_t c1c2 ;
  for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
    corrr += creal( *diq ) * creal( *s ) - cimag( *diq ) * cimag( *s ) ;
    corri += creal( *diq ) * cimag( *s ) + cimag( *diq ) * creal( *s ) ;
    diq++ , s++ ;
  }

  return corrr + I * corri;
}

// cross color 
static inline void
cross_color( double complex *__restrict a ,
	     const double complex *__restrict b ,
	     const double complex *__restrict c ,
	     const size_t id1 ,
	     const size_t id2 )
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
#elif NC == 2
  *a += b[ id1 + 0 ] * c[ id2 + 0 ] - b[ id2 + 0 ] * c[ id1 + 0 ] ; a++ ;
  *a += b[ id1 + 0 ] * c[ id2 + 2 ] - b[ id2 + 0 ] * c[ id1 + 2 ] ; a++ ;
  *a += b[ id1 + 2 ] * c[ id2 + 0 ] - b[ id2 + 2 ] * c[ id1 + 0 ] ; a++ ;
  *a += b[ id1 + 2 ] * c[ id2 + 2 ] - b[ id2 + 2 ] * c[ id1 + 2 ] ; a++ ;
#else
  size_t c1 , c2 ;
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
  struct spinor T3SNK[ NC ] ;

#if NC == 3
  // initialise temporary storage
  spinor_zero_site( &T3SNK[ 0 ] ) ;
  spinor_zero_site( &T3SNK[ 1 ] ) ;
  spinor_zero_site( &T3SNK[ 2 ] ) ;
  // Sink cross color and trace, this leaves only one set of dirac indices
  size_t i , j , d ;
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
#elif NC == 2
  fprintf( stderr , "[CROSS COLOR TRACE] NC = %d not supported\n" , NC ) ;
  fprintf( stderr , "please use diquark code instead of baryons\n" , NC ) ;
  exit(1) ;
#else
  fprintf( stderr , "[CROSS COLOR TRACE] NC = %d not supported\n" , NC ) ;
  exit(1) ;
#endif
  return ;
}

#endif
