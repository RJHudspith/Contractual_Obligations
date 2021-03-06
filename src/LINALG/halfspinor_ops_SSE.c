/**
   @file halfspinor_ops_SSE.c
   @brief vectorised halfspinor operations
 */
#include "common.h"

#include "matrix_ops.h"

#ifdef HAVE_IMMINTRIN_H

#ifdef __AVX__
  #include "AVX_OPS.h"
#endif

#if NC==3
#define inline_su3(pA,pB,pC) \
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 0 ) , *( pC + 0 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 1 ) , *( pC + 3 ) ) , \
				SSE2_MUL( *( pB + 2 ) , *( pC + 6 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 0 ) , *( pC + 1 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 1 ) , *( pC + 4 ) ) ,	\
				SSE2_MUL( *( pB + 2 ) , *( pC + 7 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 0 ) , *( pC + 2 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 1 ) , *( pC + 5 ) ) ,	\
				SSE2_MUL( *( pB + 2 ) , *( pC + 8 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 3 ) , *( pC + 0 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 4 ) , *( pC + 3 ) ) ,	\
				SSE2_MUL( *( pB + 5 ) , *( pC + 6 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 3 ) , *( pC + 1 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 4 ) , *( pC + 4 ) ) ,	\
				SSE2_MUL( *( pB + 5 ) , *( pC + 7 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 3 ) , *( pC + 2 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 4 ) , *( pC + 5 ) ) ,	\
				SSE2_MUL( *( pB + 5 ) , *( pC + 8 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 6 ) , *( pC + 0 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 7 ) , *( pC + 3 ) ) , \
				SSE2_MUL( *( pB + 8 ) , *( pC + 6 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 6 ) , *( pC + 1 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 7 ) , *( pC + 4 ) ) ,	\
				SSE2_MUL( *( pB + 8 ) , *( pC + 7 ) ) ) ) ; \
  pA++ ;								\
  *pA = _mm_add_pd( SSE2_MUL( *( pB + 6 ) , *( pC + 2 ) ) ,		\
		    _mm_add_pd( SSE2_MUL( *( pB + 7 ) , *( pC + 5 ) ) ,	\
				SSE2_MUL( *( pB + 8 ) , *( pC + 8 ) ) ) ) ; \
  pA++ ;								
#endif

// atomically add halfspinors a += b
void
add_halfspinor( struct halfspinor *a ,
		const struct halfspinor b )
{
  __m128d *pA = (__m128d*)a -> D ;
  const __m128d *pB = (const __m128d*)b.D ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *pA = _mm_add_pd( *pA , *pB ) ; pA++ ; pB++ ;
    *pA = _mm_add_pd( *pA , *pB ) ; pA++ ; pB++ ;
    *pA = _mm_add_pd( *pA , *pB ) ; pA++ ; pB++ ;
    *pA = _mm_add_pd( *pA , *pB ) ; pA++ ; pB++ ;
  }
  return ;
}

void
Fmunu_halfspinor( struct halfspinor *a ,
		  const double complex *b ,
		  const struct halfspinor c )
{
  multab( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multab( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multab( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multab( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
  return ;
}

// this thing gets called so often I really spent some time optimising it:
// the function calls _mm256d_addsub_pd a lot and _mmm256_setr_pd
// the addsub is probably unavoidable and I elimiated a lot of the setrs
// with avx2 we can actually do a better job as there is an instruction that
// allows for copies in different lanes
void
colormatrix_halfspinor( __m128d *pA ,
		        const __m128d *pB ,
		        const __m128d *pC )
{
#if (defined __AVX__) && (NC==3) && (ND==4)
  double *a = (double*)pA ;
  __m256d B[ NCNC ] ;
  const double *c = (const double*)pC ;
  // broadcast is an instruction that copies pB into upper and
  // lower parts of the 256d register, so packs (re,im,re,im) into
  // the register
  B[0] = _mm256_broadcast_pd( pB+0 ) ;
  B[1] = _mm256_broadcast_pd( pB+1 ) ;
  B[2] = _mm256_broadcast_pd( pB+2 ) ;
  B[3] = _mm256_broadcast_pd( pB+3 ) ;
  B[4] = _mm256_broadcast_pd( pB+4 ) ;
  B[5] = _mm256_broadcast_pd( pB+5 ) ;
  B[6] = _mm256_broadcast_pd( pB+6 ) ;
  B[7] = _mm256_broadcast_pd( pB+7 ) ;
  B[8] = _mm256_broadcast_pd( pB+8 ) ;
  
  // set these constants, are really just a column vector
  register __m256d c1 , c2 , c3 ;
  c1 = _mm256_loadu_pd( c ) ;
  c2 = _mm256_loadu_pd( c+6 ) ;
  c3 = _mm256_loadu_pd( c+12 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[0] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[1] , c2 ) ,
						      AVX_MUL( B[2] , c3 ) ) ) ) ;
  a += 4 ;

  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+2 ) , pC[0] , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+5 ) , pC[3] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+8 ) , pC[6] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[3] , B[0] , 3 ) , c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[4] , B[1] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[5] , B[2] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+2 ) ;
  c2 = _mm256_loadu_pd( c+8 ) ;
  c3 = _mm256_loadu_pd( c+14 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[3] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[4] , c2 ) ,
						      AVX_MUL( B[5] , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c ) ;
  c2 = _mm256_loadu_pd( c+6 ) ;
  c3 = _mm256_loadu_pd( c+12 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[6] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[7] , c2 ) ,
						      AVX_MUL( B[8] , c3 ) ) ) ) ;
  a += 4 ;
  
  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+2 ) , pC[9]  , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+5 ) , pC[12] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+8 ) , pC[15] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[0] , B[6] , 3 ) , c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[1] , B[7] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[2] , B[8] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+20 ) ;
  c2 = _mm256_loadu_pd( c+26 ) ;
  c3 = _mm256_loadu_pd( c+32 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[0] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[1] , c2 ) ,
						      AVX_MUL( B[2] , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+18 ) ;
  c2 = _mm256_loadu_pd( c+24 ) ;
  c3 = _mm256_loadu_pd( c+30 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[3] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[4] , c2 ) ,
						      AVX_MUL( B[5] , c3 ) ) ) ) ;
  a += 4 ;

  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+11 ) , pC[9]  , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+14 ) , pC[12] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+17 ) , pC[15] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[6] , B[3] , 3 ), c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[7] , B[4] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[8] , B[5] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;
    
  c1 = _mm256_loadu_pd( c+20 ) ;
  c2 = _mm256_loadu_pd( c+26 ) ;
  c3 = _mm256_loadu_pd( c+32 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[6] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[7] , c2 ) ,
						      AVX_MUL( B[8] , c3 ) ) ) ) ;
  a += 4 ;

  // increment C pointer to dirac indices 2 and 3
  c += 36 ; pC += 18 ;

  c1 = _mm256_loadu_pd( c ) ;
  c2 = _mm256_loadu_pd( c+6 ) ;
  c3 = _mm256_loadu_pd( c+12 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[0] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[1] , c2 ) ,
						      AVX_MUL( B[2] , c3 ) ) ) ) ;
  a += 4 ;

  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+2 ) , pC[0] , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+5 ) , pC[3] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+8 ) , pC[6] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[3] , B[0] , 3 ) , c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[4] , B[1] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[5] , B[2] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+2 ) ;
  c2 = _mm256_loadu_pd( c+8 ) ;
  c3 = _mm256_loadu_pd( c+14 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[3] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[4] , c2 ) ,
						      AVX_MUL( B[5] , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c ) ;
  c2 = _mm256_loadu_pd( c+6 ) ;
  c3 = _mm256_loadu_pd( c+12 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[6] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[7] , c2 ) ,
						      AVX_MUL( B[8] , c3 ) ) ) ) ;
  a += 4 ;
  
  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+2 ) , pC[9] , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+5 ) , pC[12] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+8 ) , pC[15] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[0] , B[6] , 3 ) , c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[1] , B[7] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[2] , B[8] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+20 ) ;
  c2 = _mm256_loadu_pd( c+26 ) ;
  c3 = _mm256_loadu_pd( c+32 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[0] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[1] , c2 ) ,
						      AVX_MUL( B[2] , c3 ) ) ) ) ;
  a += 4 ;

  c1 = _mm256_loadu_pd( c+18 ) ;
  c2 = _mm256_loadu_pd( c+24 ) ;
  c3 = _mm256_loadu_pd( c+30 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[3] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[4] , c2 ) ,
						      AVX_MUL( B[5] , c3 ) ) ) ) ;
  a += 4 ;

  c1 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+11 ) , pC[9] , 1 ) ;
  c2 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+14 ) , pC[12] , 1 ) ;
  c3 =  _mm256_insertf128_pd( _mm256_broadcast_pd( pC+17 ) , pC[15] , 1 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[6] , B[3] , 3 ) , c1 ) ,
				       _mm256_add_pd( AVX_MUL( _mm256_blend_pd( B[7] , B[4] , 3 ) , c2 ) ,
						      AVX_MUL( _mm256_blend_pd( B[8] , B[5] , 3 ) , c3 ) ) ) ) ;
  a += 4 ;
    
  c1 = _mm256_loadu_pd( c+20 ) ;
  c2 = _mm256_loadu_pd( c+26 ) ;
  c3 = _mm256_loadu_pd( c+32 ) ;
  _mm256_storeu_pd( a , _mm256_add_pd( AVX_MUL( B[6] , c1 ) ,
				       _mm256_add_pd( AVX_MUL( B[7] , c2 ) ,
						      AVX_MUL( B[8] , c3 ) ) ) ) ;
#else
  #if NC==3
  inline_su3( pA , pB , pC ) ; pC += NCNC ;
  inline_su3( pA , pB , pC ) ; pC += NCNC ;
  inline_su3( pA , pB , pC ) ; pC += NCNC ;
  inline_su3( pA , pB , pC ) ; 
  #else
  multab( pA , pB , pC ) ; pA += NCNC ; pC += NCNC ;
  multab( pA , pB , pC ) ; pA += NCNC ; pC += NCNC ;
  multab( pA , pB , pC ) ; pA += NCNC ; pC += NCNC ;
  multab( pA , pB , pC ) ; pA += NCNC ; pC += NCNC ;
  #endif
#endif
  return ;
}

void
colormatrixdag_halfspinor( struct halfspinor *a ,
			   const double complex b[ NCNC ] ,
			   const struct halfspinor c )
{
  multabdag( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multabdag( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multabdag( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multabdag( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
  return ;
}

void
halfspinor_Saxpy( struct halfspinor *H ,
		  const struct halfspinor S ,
		  const double fac )
{
  __m128d *pH = (__m128d*)H -> D ;
  const __m128d *pS = (const __m128d*)S.D ;
  register const __m128d f = _mm_set_pd( fac , fac ) ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *pH = SSE2_FMA( f , *pS , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , *pS , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , *pS , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , *pS , *pH ) ; pS++ ; pH++ ;
  }
  return ;
}

// version multiplied by I
void
halfspinor_iSaxpy( struct halfspinor *H ,
		   const struct halfspinor S ,
		   const double fac )
{
  __m128d *pH = (__m128d*)H -> D ;
  const __m128d *pS = (const __m128d*)S.D ;
  register const __m128d f = _mm_set_pd( fac , fac ) ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *pH = SSE2_FMA( f , SSE2_iMUL( *pS ) , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , SSE2_iMUL( *pS ) , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , SSE2_iMUL( *pS ) , *pH ) ; pS++ ; pH++ ;
    *pH = SSE2_FMA( f , SSE2_iMUL( *pS ) , *pH ) ; pS++ ; pH++ ;
  }
  return ;
}

// does H^i = H^i + \sigma^ij S^j
// where sigma_map is an index map for the sigma product
// and imap is the map in Z_4 of the elements.
// For example :
// sigma_map[ NS ] = { 2 , 3 , 0 , 1 }
// as sigma_x is { 0 , 1 , 1 , 0 }
// and imap would be { 0 , 0 , 0 , 0 } as this maps to 1
void
halfspinor_sigma_Saxpy( struct halfspinor *H ,
			const struct halfspinor S ,
			const uint8_t sigma_map[ NS ] ,
			const uint8_t imap[ NS ] )
{
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    switch( imap[d]%NS ) {
    case 0 :
      colormatrix_Saxpy( H -> D[d] , S.D[ sigma_map[d] ] , 1. ) ;
      break ;
    case 1 :
      colormatrix_iSaxpy( H -> D[d] , S.D[ sigma_map[d] ] , 1. ) ;
      break ;
    case 2 :
      colormatrix_Saxpy( H -> D[d] , S.D[ sigma_map[d] ] , -1. ) ;
      break ;
    case 3 :
      colormatrix_iSaxpy( H -> D[d] , S.D[ sigma_map[d] ] , -1. ) ;
      break ;
    }
  }
  return ;
}

// multiplies two halfspinors
// a^{alpha,beta}_{a,b} = b^{alpha,kappa}_{a,c} c^{kappa,beta}_{c,b}
void
halfspinor_multiply( struct halfspinor *a ,
		     const struct halfspinor b ,
		     const struct halfspinor c )
{
  __m128d temp[ NCNC ] ;
  // first element
  multab( (void*)a -> D[0] , (void*)b.D[0] , (void*)c.D[0] ) ;
  multab( temp , (void*)b.D[1] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[0] , temp ) ;
  // second
  multab( (void*)a -> D[1] , (void*)b.D[0] , (void*)c.D[1] ) ;
  multab( temp , (void*)b.D[1] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[1] , temp ) ;
  // third
  multab( (void*)a -> D[2] , (void*)b.D[2] , (void*)c.D[0] ) ;
  multab( temp , (void*)b.D[3] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[2] , temp ) ;
  // fourth
  multab( (void*)a -> D[3] , (void*)b.D[2] , (void*)c.D[1] ) ;
  multab( temp , (void*)b.D[3] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[3] , temp ) ;
  return ;
}

void
sigmaB_halfspinor( struct halfspinor *S1 ,
		   const struct field Fmunu ,
		   const struct halfspinor S )
{
  struct halfspinor t1 ;
  
  // initialise sigma.B into t1
  const __m128d *pF0 = (const __m128d*)Fmunu.O[0] ;
  const __m128d *pF1 = (const __m128d*)Fmunu.O[1] ;
  const __m128d *pF2 = (const __m128d*)Fmunu.O[2] ;
  __m128d *pt0 = (__m128d*)t1.D[0] ;
#if NC == 3
  // t1.D[0]
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  *pt0 = *pF2 ; pt0++ ; pF2++ ;
  // t1.D[1]
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  // t1.D[2]
  pF0 = (const __m128d*)Fmunu.O[0] ;
  pF1 = (const __m128d*)Fmunu.O[1] ;
  pF2 = (const __m128d*)Fmunu.O[2] ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  *pt0 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ; pt0++ ; pF0++ ; pF1++ ;
  // t1.D[3]
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
  *pt0 = SSE_FLIP( *( pF2 ) ) ; pt0++ ; pF2++ ;
#else
  __m128d *pt1 = (__m128d*)t1.D[1] ;
  __m128d *pt2 = (__m128d*)t1.D[2] ;
  __m128d *pt3 = (__m128d*)t1.D[3] ;
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    *pt0 = *( pF2 ) ;
    *pt1 = _mm_sub_pd( *( pF0 ) , SSE2_iMUL( *( pF1) ) ) ;
    *pt2 = _mm_add_pd( *( pF0 ) , SSE2_iMUL( *( pF1 ) ) ) ;
    *pt3 = SSE_FLIP( *( pF2 ) ) ;
    pF0++ ; pF1++ ; pF2++ ;
    pt0++ ; pt1++ ; pt2++ ; pt3++ ;
  }
#endif

  halfspinor_multiply( S1 , t1 , S ) ;
  return ;
}

// set s-matrix to zero
void
zero_halfspinor( struct halfspinor *S )
{
  __m128d *pS = (__m128d*)S -> D ;
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    *pS = _mm_setzero_pd() ; pS++ ;
    *pS = _mm_setzero_pd() ; pS++ ;
    *pS = _mm_setzero_pd() ; pS++ ;
    *pS = _mm_setzero_pd() ; pS++ ;
  }
  return ;
}

#endif
