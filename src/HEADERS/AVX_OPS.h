/**
   @file AVX_OPS.h
   @brief avx macro functions works on two double complex numbers at once
   i.e. (re,im,re,im)
 */
#ifndef AVX_OPS_H
#define AVX_OPS_H

#if !(defined __ICC)
#define SSE_FLIP(a) ( -a )
#else
#define SSE_FLIP(a) ( _mm256_xor_pd( a , _mm256_set1_pd( -0.0 ) ) )
#endif

#define AVX_MUL(a,b) ( _mm256_addsub_pd( _mm256_mul_pd( _mm256_movedup_pd( a ) , b ) , \
					 _mm256_mul_pd( _mm256_unpackhi_pd( a , a ) , \
							 _mm256_shuffle_pd( b , b , 0x5 ) ) ) )

#define AVX_MULCONJ(a,b) ( _mm256_add_pd( _mm256_mul_pd( _mm256_movedup_pd( a ) , b ) , \
					  _mm256_mul_pd( _mm256_unpackhi_pd( a , SSE_FLIP(a) ) , \
							 _mm256_shuffle_pd( b , b , 0x5 ) ) ) )

#define AVX_MUL_CONJ(a,b) ( _mm256_addsub_pd( _mm256_mul_pd( _mm256_unpackhi_pd( a , a ) , \
							     _mm256_shuffle_pd( b , b , 0x5 ) ) , \
					      _mm256_mul_pd( _mm256_movedup_pd( a ) , SSE_FLIP( b ) ) ) )

#define AVX_MUL_CONJCONJ(a,b) ( _mm256_hsub_pd( _mm256_mul_pd( a , b ) , \
						_mm256_mul_pd( a , _mm256_shuffle_pd( SSE_FLIP(b) , b , 0x5 ) ) ) )

// multiply by I
#define AVX_iMUL(a) ( _mm256_shuffle_pd( SSE_FLIP(a) , a , 0x5 ) )

// multiply by -I
#define AVX_miMUL(a) ( _mm256_shuffle_pd( a , SSE_FLIP(a) , 0x5 ) )

// complex conjugate
#define AVX_CONJ(a) ( _mm256_blend_pd( a , SSE_FLIP(a) , 0xa ) )

#endif