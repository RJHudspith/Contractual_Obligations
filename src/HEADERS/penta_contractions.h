/**
   @file penta_contractions.h
   @brief prototype declarations for pentaquark contractions
 */
#ifndef PENTA_CONTRACTIONS_H
#define PENTA_CONTRACTIONS_H

// alignments for the spinmatrix multiplies
//#if (defined __AVX512__)
//#define SPINT_ALIGNMENT (64)
//#elif (defined __AVX__)
#if (defined __AVX__)
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

#if (defined __AVX512__)
  #include "AVX512_OPS.h"
#elif (defined __AVX__)
  #include "AVX_OPS.h"
#endif

static inline void
spinmatrix_multiply_T_avx( void *a ,
			   const void *b ,
			   const void *c )
{
  /*
#if (defined __AVX512__)
  __m512d *A = (__m512d*)a ;
  const double *B = (const double*)b ;
  const double *C = (const double*)c ;

  register __m512d t1 = _mm512_setr_pd( B[0] , B[1] , B[0] , B[1] ,
					B[0] , B[1] , B[0] , B[1] ) ;
  register __m512d t2 = _mm512_setr_pd( B[2] , B[3] , B[2] , B[3] ,
					B[2] , B[3] , B[2] , B[3] ) ;
  register __m512d t3 = _mm512_setr_pd( B[4] , B[5] , B[4] , B[5] ,
					B[4] , B[5] , B[4] , B[5] ) ;
  register __m512d t4 = _mm512_setr_pd( B[6] , B[7] , B[6] , B[7] ,
					B[6] , B[7] , B[6] , B[7] ) ;

  register const __m256d t5 = _mm512_setr_pd( C[0] , C[1] , C[8] , C[9] ,
					      C[16] , C[17] , C[24] , C[25] ) ;
  register const __m512d t6 = _mm512_setr_pd( C[2] , C[3] , C[10] , C[11] ,
					      C[18] , C[19] , C[26] , C[27] ) ;
  register const __m512d t7 = _mm512_setr_pd( C[4] , C[5] , C[12] , C[13] ,
					      C[20] , C[21] , C[28] , C[29] ) ;
  register const __m512d t8 = _mm512_setr_pd( C[6] , C[7] , C[14] , C[15] ,
					      C[22] , C[23] , C[30] , C[31] ) ;
  A[0] = _mm512_add_pd( _mm512_add_pd( AVX512_MUL( t1 , t5 ) , AVX512_MUL( t2 , t6 ) ) ,
			_mm512_add_pd( AVX512_MUL( t3 , t7 ) , AVX512_MUL( t4 , t8 ) ) ) ;
  // second row
  t1 = _mm512_setr_pd( B[8] , B[9] , B[8] , B[9] , B[8] , B[9] , B[8] , B[9] ) ;
  t2 = _mm512_setr_pd( B[10] , B[11] , B[10] , B[11] , B[10] , B[11] , B[10] , B[11] ) ;
  t3 = _mm512_setr_pd( B[12] , B[13] , B[12] , B[13] , B[12] , B[13] , B[12] , B[13] ) ;
  t4 = _mm512_setr_pd( B[14] , B[15] , B[14] , B[15] , B[14] , B[15] , B[14] , B[15] ) ;
  A[1] = _mm512_add_pd( _mm512_add_pd( AVX512_MUL( t1 , t5 ) , AVX512_MUL( t2 , t6 ) ) ,
			_mm512_add_pd( AVX512_MUL( t3 , t7 ) , AVX512_MUL( t4 , t8 ) ) ) ;
  // third
  t1 = _mm512_setr_pd( B[16] , B[17] , B[16] , B[17] , B[16] , B[17] , B[16] , B[17] ) ;
  t2 = _mm512_setr_pd( B[18] , B[19] , B[18] , B[19] , B[18] , B[19] , B[18] , B[19] ) ;
  t3 = _mm512_setr_pd( B[20] , B[21] , B[20] , B[21] , B[20] , B[21] , B[20] , B[21] ) ;
  t4 = _mm512_setr_pd( B[22] , B[23] , B[22] , B[23] , B[22] , B[23] , B[22] , B[23] ) ;
  A[2] = _mm512_add_pd( _mm512_add_pd( AVX512_MUL( t1 , t5 ) , AVX512_MUL( t2 , t6 ) ) ,
			_mm512_add_pd( AVX512_MUL( t3 , t7 ) , AVX512_MUL( t4 , t8 ) ) ) ;
  // fourth
  t1 = _mm512_setr_pd( B[24] , B[25] , B[24] , B[25] , B[24] , B[25] , B[24] , B[25] ) ;
  t2 = _mm512_setr_pd( B[26] , B[27] , B[26] , B[27] , B[26] , B[27] , B[26] , B[27] ) ;
  t3 = _mm512_setr_pd( B[28] , B[29] , B[28] , B[29] , B[28] , B[29] , B[28] , B[29] ) ;
  t4 = _mm512_setr_pd( B[30] , B[31] , B[30] , B[31] , B[30] , B[31] , B[30] , B[31] ) ;
  A[3] = _mm512_add_pd( _mm512_add_pd( AVX512_MUL( t1 , t5 ) , AVX512_MUL( t2 , t6 ) ) ,
			_mm512_add_pd( AVX512_MUL( t3 , t7 ) , AVX512_MUL( t4 , t8 ) ) ) ;
#elif (defined __AVX__)
  */
#if (defined __AVX__)
  __m256d *A = (__m256d*)a ;
  const double *B = (const double*)b ;
  const double *C = (const double*)c ;

  // initially set these up into the y registers hopefully
  register __m256d t1 = _mm256_setr_pd( B[0] , B[1] , B[0] , B[1] ) ;
  register __m256d t2 = _mm256_setr_pd( B[2] , B[3] , B[2] , B[3] ) ;
  register __m256d t3 = _mm256_setr_pd( B[4] , B[5] , B[4] , B[5] ) ;
  register __m256d t4 = _mm256_setr_pd( B[6] , B[7] , B[6] , B[7] ) ;

  // first
  register const __m256d t5 = _mm256_setr_pd( C[0] , C[1] , C[8] , C[9] ) ;
  register const __m256d t6 = _mm256_setr_pd( C[2] , C[3] , C[10] , C[11] ) ;
  register const __m256d t7 = _mm256_setr_pd( C[4] , C[5] , C[12] , C[13] ) ;
  register const __m256d t8 = _mm256_setr_pd( C[6] , C[7] , C[14] , C[15] ) ;
  register const __m256d t9 = _mm256_setr_pd( C[16] , C[17] , C[24] , C[25] ) ;
  register const __m256d t10 = _mm256_setr_pd( C[18] , C[19] , C[26] , C[27] ) ;
  register const __m256d t11 = _mm256_setr_pd( C[20] , C[21] , C[28] , C[29] ) ;
  register const __m256d t12 = _mm256_setr_pd( C[22] , C[23] , C[30] , C[31] ) ;
  
  A[0] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[1] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) ,  AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) , AVX_MUL( t4 , t12 ) ) ) ;
  // second
  t1 = _mm256_setr_pd( B[8] , B[9] , B[8] , B[9] ) ;
  t2 = _mm256_setr_pd( B[10] , B[11] , B[10] , B[11] ) ;
  t3 = _mm256_setr_pd( B[12] , B[13] , B[12] , B[13] ) ;
  t4 = _mm256_setr_pd( B[14] , B[15] , B[14] , B[15] ) ;
  A[2] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[3] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) , AVX_MUL( t2 , t10 ) ) ,
			_mm256_add_pd( AVX_MUL( t3 , t11 ) ,  AVX_MUL( t4 , t12 ) ) ) ;
  // third
  t1 = _mm256_setr_pd( B[16] , B[17] , B[16] , B[17] ) ;
  t2 = _mm256_setr_pd( B[18] , B[19] , B[18] , B[19] ) ;
  t3 = _mm256_setr_pd( B[20] , B[21] , B[20] , B[21] ) ;
  t4 = _mm256_setr_pd( B[22] , B[23] , B[22] , B[23] ) ;
  A[4] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t5 ) , AVX_MUL( t2 , t6 ) ) ,
		      _mm256_add_pd( AVX_MUL( t3 , t7 ) , AVX_MUL( t4 , t8 ) ) ) ;
  A[5] = _mm256_add_pd( _mm256_add_pd( AVX_MUL( t1 , t9 ) , AVX_MUL( t2 , t10 ) ) ,
		      _mm256_add_pd( AVX_MUL( t3 , t11 ) ,  AVX_MUL( t4 , t12 ) ) ) ;
  // fourth
  t1 = _mm256_setr_pd( B[24] , B[25] , B[24] , B[25] ) ;
  t2 = _mm256_setr_pd( B[26] , B[27] , B[26] , B[27] ) ;
  t3 = _mm256_setr_pd( B[28] , B[29] , B[28] , B[29] ) ;
  t4 = _mm256_setr_pd( B[30] , B[31] , B[30] , B[31] ) ;
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
