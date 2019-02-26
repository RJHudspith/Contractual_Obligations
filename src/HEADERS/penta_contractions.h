/**
   @file penta_contractions.h
   @brief prototype declarations for pentaquark contractions
 */
#ifndef PENTA_CONTRACTIONS_H
#define PENTA_CONTRACTIONS_H

// alignments for the spinmatrix multiplies
#if (defined __AVX__) && (defined HAVE_IMMINTRIN_H)
#define SPINT_ALIGNMENT (32)
#else
#define SPINT_ALIGNMENT (16)
#endif

/**
   @fn static inline size_t idx( const size_t b , const size_t bp , const size_t c , const size_t cp , const size_t g , const size_t gp , const size_t h , const size_t hp )
   @brief return a linearised index from a tensor with 8 indices
 */
static inline size_t
idx( const size_t b , const size_t bp ,
     const size_t c , const size_t cp ,
     const size_t g , const size_t gp ,
     const size_t h , const size_t hp )
{
  return b + NC * ( bp + NC * ( c + NC * ( cp + NC * ( g + NC * ( gp + NC * ( h + NC * hp ) ) ) ) ) ) ;
}

/**
   @fn static inline size_t idx( const size_t b , const size_t bp , const size_t c , const size_t cp , const size_t g , const size_t gp , const size_t h , const size_t hp )
   @brief return a linearised index from a tensor with 6 indices
 */
static inline size_t
idx2( const size_t b , const size_t bp ,
      const size_t c , const size_t cp ,
      const size_t g , const size_t gp )
{
  return b + NC * ( bp + NC * ( c + NC * ( cp + NC * ( g + NC * gp ) ) ) ) ;
}

#if (defined __AVX__) && (HAVE_IMMINTRIN_H)
  #include "AVX_OPS.h"
#else
  #include "spinmatrix_ops.h"
#endif

static inline void
spinmatrix_multiply_T_avx( void *a ,
			   const void *b ,
			   const void *c )
{
#if (defined __AVX__) && (defined HAVE_IMMINTRIN_H)
  __m256d *A = (__m256d*)a ;
  const __m128d *B = (const __m128d*)b ;
  const __m128d *C = (const __m128d*)c ;

  // first bunch of variables
  register const __m256d t5 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+0 )   , C[4]  , 1 ) ;
  register const __m256d t6 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+1 )   , C[5]  , 1 ) ;
  register const __m256d t7 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+2 )   , C[6]  , 1 ) ;
  register const __m256d t8 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+3 )   , C[7]  , 1 ) ;
  register const __m256d t9 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+8 )   , C[12] , 1 ) ;
  register const __m256d t10 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+9 )  , C[13] , 1 ) ;
  register const __m256d t11 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+10 ) , C[14] , 1 ) ;
  register const __m256d t12 = _mm256_insertf128_pd( _mm256_broadcast_pd( C+11 ) , C[15] , 1 ) ;

  // initially set these up into the y registers hopefully and do the multiply
  register __m256d t1 = _mm256_broadcast_pd( B+0 ) ;
  register __m256d t2 = _mm256_broadcast_pd( B+1 ) ;
  register __m256d t3 = _mm256_broadcast_pd( B+2 ) ;
  register __m256d t4 = _mm256_broadcast_pd( B+3 ) ;
  A[0] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[1] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) ,  AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) , AVX_MUL( t4 , t12 ) ) ) ;
  // second
  t1 = _mm256_broadcast_pd( B+4 ) ;
  t2 = _mm256_broadcast_pd( B+5 ) ;
  t3 = _mm256_broadcast_pd( B+6 ) ;
  t4 = _mm256_broadcast_pd( B+7 ) ;
  A[2] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[3] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) , AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) ,  AVX_MUL( t4 , t12 ) ) ) ;
  // third
  t1 = _mm256_broadcast_pd( B+8 ) ;
  t2 = _mm256_broadcast_pd( B+9 ) ;
  t3 = _mm256_broadcast_pd( B+10 ) ;
  t4 = _mm256_broadcast_pd( B+11 ) ;
  A[4] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[5] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) , AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) ,  AVX_MUL( t4 , t12 ) ) ) ;
  // fourth
  t1 = _mm256_broadcast_pd( B+12 ) ;
  t2 = _mm256_broadcast_pd( B+13 ) ;
  t3 = _mm256_broadcast_pd( B+14 ) ;
  t4 = _mm256_broadcast_pd( B+15 ) ;
  A[6] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[7] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) , AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) , AVX_MUL( t4 , t12 ) ) ) ;
#else
  spinmatrix_multiply_T( a , b , c ) ;
#endif
  return ;
}

/**
   @fn void spinmatrix_multiply_T_avx( void *a , const void *b , const void *c )
   @brief computes a = b*c^T with the fastest instruction set possible
 */
void
spinmatrix_multiply_T_avx( void *a ,
			   const void *b ,
			   const void *c ) ;

/**
   @fn int pentas( double complex *result , const struct spinor L , const struct spinor S , const struct spinor bwdH , const struct gamma *GAMMAS )
   @brief pentaquark contraction code for a udusb pentaquark
   @param L :: light quark propagator assumes u-d degeneracy
   @param S :: strange quark propagator
   @param bwdH :: backward-propagating heavy quark
   @param GAMMAS :: gamma matrices
   @param loc :: indexing
 */
int
pentas( double complex *result ,
	double complex **F ,
	const struct spinor L , 
	const struct spinor S ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS ,
	const uint8_t **loc ) ;

#endif
