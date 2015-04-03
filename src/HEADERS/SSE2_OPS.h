/**
   @file SSE2_OPS.h
   @brief inline macros for SSE2 operations
 */
#ifndef SSE2_OPS_H
#define SSE2_OPS_H

#ifdef __SSE3__
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , -a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , -a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

#ifdef __SSE3__
#define SSE2_MUL(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
				       _mm_mul_pd( _mm_unpackhi_pd( a , a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MUL(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
				    _mm_mul_pd( _mm_unpackhi_pd( -a , a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

#define SSE2_MUL_CONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , -a ) , b ) , \
					 _mm_mul_pd( _mm_unpackhi_pd( a , a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )

#define SSE2_MUL_CONJCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , -a ) , b ) , \
					     _mm_mul_pd( _mm_unpackhi_pd( -a , -a ) , _mm_shuffle_pd( b , b , 1 ) ) ) )

#define SSE2_iMUL(a) ( _mm_shuffle_pd( -a , a , 1 ) )

#define SSE2_miMUL(a) ( _mm_shuffle_pd( a , -a , 1 ) )

#define SSE2_CONJ(a) ( _mm_move_sd( -a, a ) ) ;

#endif
