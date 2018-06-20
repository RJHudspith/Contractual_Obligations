/**
   @file halfspinor_ops.c
   @brief perform arithmetic operations on halfspinors
 */
#include "common.h"

#include "matrix_ops.h"

#ifdef HAVE_EMMINTRIN_H

// pA = pB.pC where pB is hermitian and pC is an ordinary 3x3 matrix
static inline void
herm_multab( __m128d *__restrict pA ,
	     const __m128d *__restrict pB ,
	     const __m128d *__restrict pC )
{
#if NC == 3 
  register const __m128d reA = _mm_movedup_pd( *( pB + 0 ) ) ;
  register const __m128d reD = _mm_movedup_pd( *( pB + 4 ) ) ;
  register const __m128d reF = SSE_FLIP( _mm_add_pd( reA , reD ) ) ;
  // first row
  *pA = _mm_add_pd( _mm_mul_pd( reA , *( pC + 0 ) ) , 
		    SSE2_MUL( *( pB + 1 ) , *( pC + 3 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 2 ) , *( pC + 6 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( _mm_mul_pd( reA , *( pC + 1 ) ) , 
		    SSE2_MUL( *( pB + 1 ) , *( pC + 4 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 2 ) , *( pC + 7 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( _mm_mul_pd( reA , *( pC + 2 ) ) , 
		    SSE2_MUL( *( pB + 1 ) , *( pC + 5 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 2 ) , *( pC + 8 ) ) ) ; pA++ ;
  // sepCond row
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 1 ) , *( pC + 0 ) ) , 
		    _mm_mul_pd( reD , *( pC + 3 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 5 ) , *( pC + 6 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 1 ) , *( pC + 1 ) ) , 
		    _mm_mul_pd( reD , *( pC + 4 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 5 ) , *( pC + 7 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 1 ) , *( pC + 2 ) ) , 
		    _mm_mul_pd( reD , *( pC + 5 ) ) ) ;
  *pA = _mm_add_pd( *pA , SSE2_MUL( *( pB + 5 ) , *( pC + 8 ) ) ) ; pA++ ;
  // third row
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 2 ) , *( pC + 0 ) ) , 
		    SSE2_MULCONJ( *( pB + 5 ) , *( pC + 3 ) ) ) ;
  *pA = _mm_add_pd( *pA , _mm_mul_pd( reF , *( pC + 6 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 2 ) , *( pC + 1 ) ) , 
		    SSE2_MULCONJ( *( pB + 5 ) , *( pC + 4 ) ) ) ;
  *pA = _mm_add_pd( *pA , _mm_mul_pd( reF , *( pC + 7 ) ) ) ; pA++ ;
  *pA = _mm_add_pd( SSE2_MULCONJ( *( pB + 2 ) , *( pC + 2 ) ) , 
		    SSE2_MULCONJ( *( pB + 5 ) , *( pC + 5 ) ) ) ;
  *pA = _mm_add_pd( *pA , _mm_mul_pd( reF , *( pC + 8 ) ) ) ;
#else
  multab( pA , pB , pC ) ;
#endif
}
#else
static inline void
herm_multab( double complex *a ,
	     const double complex *b ,
	     const double complex *c )
{
  multab( a , b , c ) ;
  return ;
}
#endif

// atomically add halfspinors a += b
void
add_halfspinor( struct halfspinor *a ,
		const struct halfspinor b )
{
  add_mat( (void*)a -> D[0] , (void*)b.D[0] ) ;
  add_mat( (void*)a -> D[1] , (void*)b.D[1] ) ;
  add_mat( (void*)a -> D[2] , (void*)b.D[2] ) ;
  add_mat( (void*)a -> D[3] , (void*)b.D[3] ) ;
  return ;
}

void
Fmunu_halfspinor( struct halfspinor *a ,
		  const double complex *b ,
		  const struct halfspinor c )
{
  herm_multab( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  herm_multab( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  herm_multab( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  herm_multab( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
  return ;
}

void
colormatrix_halfspinor( struct halfspinor *a ,
			const double complex b[ NCNC ] ,
			const struct halfspinor c )
{
  multab( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multab( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multab( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multab( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
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
  colormatrix_Saxpy( H -> D[0] , S.D[0] , fac ) ;
  colormatrix_Saxpy( H -> D[1] , S.D[1] , fac ) ;
  colormatrix_Saxpy( H -> D[2] , S.D[2] , fac ) ;
  colormatrix_Saxpy( H -> D[3] , S.D[3] , fac ) ;
  return ;
}

// version multiplied by I
void
halfspinor_iSaxpy( struct halfspinor *H ,
		   const struct halfspinor S ,
		   const double fac )
{
  colormatrix_iSaxpy( H -> D[0] , S.D[0] , fac ) ;
  colormatrix_iSaxpy( H -> D[1] , S.D[1] , fac ) ;
  colormatrix_iSaxpy( H -> D[2] , S.D[2] , fac ) ;
  colormatrix_iSaxpy( H -> D[3] , S.D[3] , fac ) ;
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

void
halfspinor_multiply( struct halfspinor *a ,
		     const struct halfspinor b ,
		     const struct halfspinor c )
{
  double complex temp[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  // first element
  multab( (void*)a -> D[0] , (void*)b.D[0] , (void*)c.D[0] ) ;
  multab( (void*)temp , (void*)b.D[1] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[0] , (void*)temp ) ;
  // second
  multab( (void*)a -> D[1] , (void*)b.D[0] , (void*)c.D[1] ) ;
  multab( (void*)temp , (void*)b.D[1] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[1] , (void*)temp ) ;
  // third
  multab( (void*)a -> D[2] , (void*)b.D[2] , (void*)c.D[0] ) ;
  multab( (void*)temp , (void*)b.D[3] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[2] , (void*)temp ) ;
  // fourth
  multab( (void*)a -> D[3] , (void*)b.D[2] , (void*)c.D[1] ) ;
  multab( (void*)temp , (void*)b.D[3] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[3] , (void*)temp ) ;
  return ;
}

// set s-matrix to zero
void
zero_halfspinor( struct halfspinor *S )
{
  zero_colormatrix( S -> D[0] ) ;
  zero_colormatrix( S -> D[1] ) ;
  zero_colormatrix( S -> D[2] ) ;
  zero_colormatrix( S -> D[3] ) ;
  return ;
}
