/**
   @file bar_contractions.c
   @brief baryon contraction codes
 */

#include "common.h"

#include "contractions.h" // gamma_mul_lr()
#include "spinor_ops.h"   // spinor_zero_site()

// This contracts the diquark with the remaining propagator
// This does the color trace Tr[ A . B ] 
static const double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const int d0 ,
		 const int d1 ,
		 const int d2 ,
		 const int d3 )
{
#if (defined HAVE_IMMINTRIN_H) && (NC == 3)
  const __m128d *d = (const __m128d*)DiQ.D[d1][d0].C ;
  const __m128d *s = (const __m128d*)S.D[d2][d3].C ;
  register __m128d sum = _mm_setzero_pd() ;
  sum = SSE2_MUL( *d , *s) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  sum = _mm_add_pd( sum , SSE2_MUL( *d , *s ) ) ; d++ ; s++ ;
  double complex res ;
  _mm_store_pd( (void*)&res , sum ) ;
  return res ;
#else
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
#endif
}

// This carries out the color cross product and traces one set of Dirac indices.
// The results forms a diquark-type object
static void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S )
{
  // temporary 3-spinor space
  struct spinor T3SNK[ 3 ] ;
  spinor_zero_site( &T3SNK[ 0 ] ) ;
  spinor_zero_site( &T3SNK[ 1 ] ) ;
  spinor_zero_site( &T3SNK[ 2 ] ) ;

  int c1, c2, i, j, d;
  // Sink cross color and trace, this leaves only one set of dirac indices
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      // Dirac sink trace
      for( d = 0 ; d < NS ; d++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    T3SNK[0].D[i][j].C[c1][c2] += S.D[i][d].C[c1][1] * DiQ -> D[j][d].C[c2][2] ; 
	    T3SNK[0].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][2] * DiQ -> D[j][d].C[c2][1] ;
	    T3SNK[1].D[i][j].C[c1][c2] += S.D[i][d].C[c1][2] * DiQ -> D[j][d].C[c2][0] ;
	    T3SNK[1].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][0] * DiQ -> D[j][d].C[c2][2] ;
	    T3SNK[2].D[i][j].C[c1][c2] += S.D[i][d].C[c1][0] * DiQ -> D[j][d].C[c2][1] ;
	    T3SNK[2].D[i][j].C[c1][c2] -= S.D[i][d].C[c1][1] * DiQ -> D[j][d].C[c2][0] ;
	  }
	}
      }
    }
  }

  // SSE2 version of the below
#ifdef HAVE_IMMINTRIN_H
  __m128d *D1 = (__m128d*)DiQ -> D ;
  const __m128d *t0 = (const __m128d*)T3SNK[0].D ;
  const __m128d *t1 = (const __m128d*)T3SNK[1].D ;
  const __m128d *t2 = (const __m128d*)T3SNK[2].D ;
  for( i = 0 ; i < NSNS ; i++ ) {
    *D1 = _mm_sub_pd( *( t0 + NC + 2 ) , *( t0 + 1 + 2 * NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + NC + 2 ) , *( t1 + 1 + 2 * NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + NC + 2 ) , *( t2 + 1 + 2 * NC ) ) ; D1++ ;

    *D1 = _mm_sub_pd( *( t0 + 2 * NC ) , *( t0 + 2 ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + 2 * NC ) , *( t1 + 2 ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + 2 * NC ) , *( t2 + 2 ) ) ; D1++ ;

    *D1 = _mm_sub_pd( *( t0 + 1 ) , *( t0 + NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t1 + 1 ) , *( t1 + NC ) ) ; D1++ ;
    *D1 = _mm_sub_pd( *( t2 + 1 ) , *( t2 + NC ) ) ; D1++ ;

    t0 += NCNC ; t1 += NCNC ; t2 += NCNC ;
  }
#else
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
#endif
  return ;
}

// baryon contraction of ( \Gamma ) ( \Gamma^{T} ) S3 ( S2 S1 )
void
baryon_contract_site( double complex **term ,
		      const struct spinor S1 , 
		      const struct spinor S2 , 
		      const struct spinor S3 , 
		      const struct gamma Cgmu ,
		      const struct gamma CgmuT )
{
  // get diquark
  struct spinor DiQ = S1 ;
  gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;

  // Cross color product and sink Dirac trace back into DiQ
  cross_color_trace( &DiQ , S2 ) ;

  // loop over open dirac indices
  int odc , OD1 , OD2 ;
  for( odc = 0 ; odc < NSNS ; odc++ ) {
    // open dirac index for source and sink
    OD1 = odc / NS ;
    OD2 = odc % NS ;
    // Contract with the final propagator and trace out the source Dirac indices
    // A polarization must still be picked for the two open Dirac indices offline
#if NS == 4
    term[0][ odc ] += baryon_contract( DiQ, S3 , 0 , 0 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 1 , 1 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 2 , 2 , OD1 , OD2 ) ;
    term[0][ odc ] += baryon_contract( DiQ, S3 , 3 , 3 , OD1 , OD2 ) ;

    term[1][ odc ] += baryon_contract( DiQ, S3 , 0 , OD1 , 0 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 1 , OD1 , 1 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 2 , OD1 , 2 , OD2 ) ;
    term[1][ odc ] += baryon_contract( DiQ, S3 , 3 , OD1 , 3 , OD2 ) ; 
#else
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      term[0][ odc ] += baryon_contract( DiQ, S3 , dirac , dirac , OD1 , OD2 ) ;
      term[1][ odc ] += baryon_contract( DiQ, S3 , dirac , OD1 , dirac , OD2 ) ;
    }
#endif
  }
  return ;
}

// baryon contraction for the omega
void
baryon_contract_omega_site( double complex **term ,
			    const struct spinor S1 , 
			    const struct spinor S2 , 
			    const struct spinor S3 , 
			    const struct gamma Cgmu ,
			    const struct gamma CgmuT )
{
  // get diquark
  struct spinor DiQ = S1 , CgS51 = S1 , CgS52 = S1 ;
  gamma_mul_lr( &DiQ , CgmuT , Cgmu ) ;
  gamma_mul_l( &CgS51 , CgmuT ) ;
  gamma_mul_r( &CgS52 , Cgmu ) ;

  // Cross color product and sink Dirac trace back into DiQ
  struct spinor DiQ51 = CgS51 , DiQ52 = CgS52 ;
  cross_color_trace( &DiQ , S2 ) ;
  cross_color_trace( &DiQ51 , CgS52 ) ;
  cross_color_trace( &DiQ52 , S2 ) ;

  // loop over open dirac indices
  int odc , OD1 , OD2 ;
  for( odc = 0 ; odc < NSNS ; odc++ ) {
    // open dirac index for source and sink
    OD1 = odc / NS ;
    OD2 = odc % NS ;
    // Contract with the final propagator and trace out the source Dirac indices
    // A polarization must still be picked for the two open Dirac indices offline
    int dirac ;
    for( dirac = 0 ; dirac < NS ; dirac++ ){
      term[0][odc] += baryon_contract( DiQ , S3 , dirac , dirac , OD1 , OD2 ) ;
      term[1][odc] += baryon_contract( DiQ , S3 , dirac , OD1 , dirac , OD2 ) ;
      term[2][odc] += baryon_contract( DiQ51 , S3 , dirac , dirac , OD1 , OD2 ) ;
      term[3][odc] += baryon_contract( DiQ52 , CgS52, dirac , OD1 , dirac , OD2 ) ;
      term[4][odc] += baryon_contract( DiQ52 , CgS52, OD1 , dirac , dirac , OD2 ) ;
      term[5][odc] += baryon_contract( DiQ51 , S3 , OD1 , dirac , dirac , OD2 ) ;
    }
  }
  return ;
}
