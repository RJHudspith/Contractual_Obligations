#ifndef AVX512_OPS_H
#define AVX512_OPS_H

// complex multiply
#define AVX512_MUL(a,b)( _mm512_fmaddsub_pd( _mm512_movedup_pd( A ) , B , \
					     _mm512_mul_pd( _mm512_unpackhi_pd( A , A ) , \
							    _mm512_permute_pd( B , _MM_SHUFFLE(1,1,1,1) ))) )

#endif
