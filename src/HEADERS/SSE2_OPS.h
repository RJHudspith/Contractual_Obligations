/**
   @file SSE2_OPS.h
   @brief inline macros for SSE2 or SSE3 operations

   I use the first element in the __m128d structure as the real part and the second
   element as the imaginary part

   If you want to check this you pack an __m128d object with a complex a
   __m128d ma = _mm_setr_pd( creal(a) , cimag(a) )

   Some of the SSE3 intrinsics have a shorter number of instructions compiled under gcc
   so I would recommend using -msse3 instead of -msse2 wherever possible

   These are macros which are expanded in the code to ensure inlining
 */
#ifndef SSE2_OPS_H
#define SSE2_OPS_H

// gcc/clang allow for + / - * with SSE types, icc does not
#ifdef __GNUC__
  #define SSE_FLIP(a) ( -a )
#else
  #define SSE_FLIP(a) ( _mm_xor_pd( a , _mm_set1_pd( -0.0 ) ) )
#endif

// performs conj(a) * b using intrinsics
#ifdef __SSE3__
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						    _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MULCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
					_mm_mul_pd( _mm_unpackhi_pd( a , SSE_FLIP(a) ) , \
						    _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * b using intrinsics
#ifdef __SSE3__
#define SSE2_MUL(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_movedup_pd( a ) , b ) , \
				       _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						   _mm_shuffle_pd( b , b , 1 ) ) ) )
#else
#define SSE2_MUL(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) , \
				    _mm_mul_pd( _mm_unpackhi_pd( SSE_FLIP(a) , a ) , \
						_mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs a * conj( b ) using intrinsics
#ifdef __SSE3__
#define SSE2_MUL_CONJ(a,b) ( _mm_addsub_pd( _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
							_mm_shuffle_pd( b , b , 1 ) ) , \
					    _mm_mul_pd( _mm_movedup_pd( a ) , SSE_FLIP( b ) ) ) )
#else
#define SSE2_MUL_CONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , SSE_FLIP(a) ) , b ) , \
					 _mm_mul_pd( _mm_unpackhi_pd( a , a ) , \
						     _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// performs conj( a ) * conj( b ) using intrinsics
#ifdef __SSE3__
#define SSE2_MUL_CONJCONJ(a,b) ( _mm_hsub_pd( _mm_mul_pd( a , b ) , \
					      _mm_mul_pd( a , _mm_shuffle_pd( SSE_FLIP(b) , b , 1 ) ) ) )
#else
#define SSE2_MUL_CONJCONJ(a,b) ( _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , SSE_FLIP(a) ) , b ) , \
					     _mm_mul_pd( _mm_unpackhi_pd( SSE_FLIP(a) , SSE_FLIP(a) ) , \
							 _mm_shuffle_pd( b , b , 1 ) ) ) )
#endif

// multiply by I
#define SSE2_iMUL(a) ( _mm_shuffle_pd( SSE_FLIP(a) , a , 1 ) )

// multiply by -I
#define SSE2_miMUL(a) ( _mm_shuffle_pd( a , SSE_FLIP(a) , 1 ) )

// complex conjugate
#define SSE2_CONJ(a) ( _mm_move_sd( SSE_FLIP(a) , a ) ) ;

#endif
