/**
   @file contractions_SSE.c
   @brief contraction codes and brute force gamma multiplies (SSE2 version)
 */

#include "common.h"

#ifdef HAVE_EMMINTRIN_H

#include "contractions.h" // so we can alphabetise
#include "matrix_ops.h"   // colortrace_prod

// conjugate transpose of dirac indices
void
adjoint_spinor( struct spinor *__restrict adj ,
		const struct spinor S )
{
  size_t d1 , d2 ;
  __m128d *res = (__m128d*)adj -> D ;
  const __m128d *s ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      s = (const __m128d*)S.D[d2][d1].C ;
      #if NC == 3
      *res = SSE2_CONJ( *( s + 0 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 3 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 6 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 1 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 4 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 7 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 2 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 5 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 8 ) ) ; res++ ;
      #elif NC == 2
      *res = SSE2_CONJ( *( s + 0 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 2 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 1 ) ) ; res++ ;
      *res = SSE2_CONJ( *( s + 3 ) ) ; res++ ;
      #else
      size_t c1 , c2 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  *res = SSE2_CONJ( s[ c1 + c2 * NC ] ) ; res++ ;
	}
      }
      #endif
    }
  }
  return ;
}

// returns spin-color trace : Tr[ A B ]
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B )
{
  size_t d1 , d2 ;
  const __m128d *a = (const __m128d*)A.D ;
  const __m128d *b ;
  register __m128d sum = _mm_setzero_pd() ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      b = (__m128d*)B.D[ d2 ][ d1 ].C ;
      sum = _mm_add_pd( sum , colortrace_prod( a , b ) ) ;
      a += NCNC ;
      
    }
  }
  double complex s ;
  _mm_store_pd( (void*)&s , sum ) ;
  return s ;
}

// computes ( G5 S G5 )^{ \dagger }
void
full_adj( struct spinor *__restrict adj ,
	  const struct spinor S ,
	  const struct gamma G5 )
{
  struct spinor tmp = S ;
  gamma_mul_lr( &tmp , G5 , G5 ) ;     // left multiply by gamma_5
  adjoint_spinor( adj , tmp ) ;  // daggers a spinor
  return ;
}

// atomic left multiply by a gamma matrix
void
gamma_mul_l( struct spinor *__restrict res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  size_t i , j ;
  // not allowed to declare in a switch
#if NC > 3
  size_t c ;
#endif
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col = GAMMA.ig[i] ;
    __m128d *tp1 = (__m128d*)tmp.D[i] ;
    const __m128d *tp2 = (__m128d*)res -> D[col] ;
    switch( GAMMA.g[i] ) {
    case 0 : // do nothing
      for( j = 0 ; j < NS ; j++ ) {
	#if NC == 3 
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = *tp2 ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 1 :     // multiply by I
      for( j = 0 ; j < NS ; j++ ) {
	#if NC == 3 
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 2 : 
      for( j = 0 ; j < NS ; j++ ) {
	#if NC == 3 
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 3 : 
       for( j = 0 ; j < NS ; j++ ) {
	#if NC == 3 
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    }
  }
  //
  *res = tmp ;
}

// multiply a spinor on the right with a gamma matrix
void
gamma_mul_r( struct spinor *__restrict res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  __m128d *tp1 ;
  const __m128d *tp2 ;
  size_t i , j ;
  // can't declare this in a switch
#if NC > 3
  size_t c ;
#endif
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const int col = GAMMA.ig[j] ;
    switch( GAMMA.g[ col ] ) {
    case 0 : // copy
      for( i = 0 ; i < NS ; i++ ) {
	tp1 = (__m128d*)tmp.D[i][j].C ;
	tp2 = (__m128d*)res -> D[i][col].C ;
	#if NC == 3
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = *tp2 ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 1 : // multiply by i
      for( i = 0 ; i < NS ; i++ ) {
	tp1 = (__m128d*)tmp.D[i][j].C ;
	tp2 = (__m128d*)res -> D[i][col].C ;
	#if NC == 3
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 2 : // multiply by -1
      for( i = 0 ; i < NS ; i++ ) {
	tp1 = (__m128d*)tmp.D[i][j].C ;
	tp2 = (__m128d*)res -> D[i][col].C ;
	#if NC == 3
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
      }
      break ;
    case 3 : // multiply by -i
      for( i = 0 ; i < NS ; i++ ) {
	tp1 = (__m128d*)tmp.D[i][j].C ;
	tp2 = (__m128d*)res -> D[i][col].C ;
	#if NC == 3
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ; 
	}
	#endif
      }
      break ;
    }
  }
  *res = tmp ;
  return ;
}

//
void
gamma_mul_lr( struct spinor *__restrict S , 
	      const struct gamma GLEFT ,
	      const struct gamma GRIGHT )
{
  struct spinor tmp = *S ; // temporary space
  __m128d *tp1 ;
  const __m128d *tp2 ;
  size_t i , j , col1 , col2 ;
  // not allowed to declare in a switch
#if NC > 3
  size_t c ;
#endif
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    col1 = GLEFT.ig[ i ] ;
    
    for( j = 0 ; j < NS ; j++ ) {
      col2 = GRIGHT.ig[ j ] ;

      // __m128d pointer cache
      tp1 = (__m128d*)tmp.D[i][j].C ;
      tp2 = (const __m128d*)S -> D[col1][col2].C ;

      switch( ( GLEFT.g[ i ] + GRIGHT.g[ col2 ] ) & 3 ) {
      case 0 :
	#if NC == 3 
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	*tp1 = *tp2 ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = *tp2 ; tp1++ ; tp2++ ;
	}
	#endif
	break ;
      case 1 :
	#if NC == 3
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_iMUL( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
	break ;
      case 2 :
	#if NC == 3 
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE_FLIP( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
	break ;
      case 3 :
	#if NC == 3
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#elif NC == 2
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	*tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	#else
	for( c = 0 ; c < NCNC ; c++ ) {
	  *tp1 = SSE2_miMUL( *tp2 ) ; tp1++ ; tp2++ ;
	}
	#endif
	break ;
      }
    }
  }
  //
  *S = tmp ;
  return ;
}

// meson contraction code computes -Tr[ GSNK ( G5 bwd G5 )^{\dagger} GSRC ( fwd ) ]
double complex 
meson_contract( const struct gamma GSNK ,
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 )
{
  // global sum gets cast to double complex
  register __m128d gsum = _mm_setzero_pd( ) ;

  const register __m128d *d1 ;
  const register __m128d *d2 = (const __m128d*)fwd.D ;

  // local sum
  register __m128d sum ;
#if NC == 3
  register __m128d sum1 , sum2 ;
#endif

  size_t i , j ;
  register uint8_t col2 , col1 , G5GSRC ;
  for( j = 0 ; j < NS ; j++ ) {

    col2 = GSRC.ig[ G5.ig[ j ] ] ; 
    
    G5GSRC = G5.g[ col2 ] + GSRC.g[ col2 ] ;

    for( i = 0 ; i < NS ; i++ ) {

      col1 = G5.ig[ GSNK.ig[ i ] ] ;
      
      // cache the second color matrix
      d1 = (const __m128d*)bwd.D[col2][col1].C ;

      // unrolled for SU(3)
      #if NC == 3
      sum = SSE2_MULCONJ( *d1 , *d2 ) ;
      sum1 = _mm_add_pd( _mm_add_pd( SSE2_MULCONJ( *(d1+1) , *(d2+1) ) ,
				     SSE2_MULCONJ( *(d1+2) , *(d2+2) ) ) ,
			 _mm_add_pd( SSE2_MULCONJ( *(d1+3) , *(d2+3) ) ,
				     SSE2_MULCONJ( *(d1+4) , *(d2+4) ) ) ) ; 
      sum2 = _mm_add_pd( _mm_add_pd( SSE2_MULCONJ( *(d1+5) , *(d2+5) ) ,
				     SSE2_MULCONJ( *(d1+6) , *(d2+6) ) ) ,
			 _mm_add_pd( SSE2_MULCONJ( *(d1+7) , *(d2+7) ) ,
				     SSE2_MULCONJ( *(d1+8) , *(d2+8) ) ) ) ;
      sum = _mm_add_pd( sum , _mm_add_pd( sum1 , sum2 ) ) ; d2 += 9 ;
      #elif NC == 2
      sum = SSE2_MULCONJ( *d1 , *d2 ) ; d1 ++ ; d2 ++ ;
      sum = _mm_add_pd( sum , SSE2_MULCONJ( *d1 , *d2 ) ) ; d1 ++ ; d2 ++ ;
      sum = _mm_add_pd( sum , SSE2_MULCONJ( *d1 , *d2 ) ) ; d1 ++ ; d2 ++ ;
      sum = _mm_add_pd( sum , SSE2_MULCONJ( *d1 , *d2 ) ) ; d1 ++ ; d2 ++ ;
      #else 
      sum = _mm_setzero_pd( ) ;
      size_t c1c2 ;
      for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
	sum = _mm_add_pd( sum , SSE2_MULCONJ( *d1 , *d2 ) ) ; d1 ++ ; d2 ++ ;
      }
      #endif

      // switch for the phases (note the implicit minus sign!)
      switch( ( GSNK.g[ i ] + G5.g[ col1 ] + G5GSRC ) & 3 ) {
      case 0 : gsum = _mm_sub_pd( gsum , sum ) ; break ;
      case 1 : gsum = _mm_add_pd( gsum , _mm_shuffle_pd( sum , SSE_FLIP( sum ) , 1 ) ) ; break ;
      case 2 : gsum = _mm_add_pd( gsum , sum ) ; break ;
      case 3 : gsum = _mm_add_pd( gsum , _mm_shuffle_pd( SSE_FLIP( sum ) , sum , 1 ) ) ; break ;
      }
      // and we are done
    }
  }
  // cast through the void
  double complex s ;
  _mm_store_pd( (void*)&s , gsum ) ;
  return s ;
}

// meson contraction code computes -Tr[ GSNK ( bwd ) GSRC ( fwd ) ]
double complex
simple_meson_contract( const struct gamma GSNK ,		
		       const struct spinor bwd , 
		       const struct gamma GSRC ,
		       const struct spinor fwd )
{
  register __m128d gsum = _mm_setzero_pd( ) ;

  const __m128d *fcache ;
  const __m128d *bcache ;
  register __m128d sum ;

  size_t i , j ;
  // loop columns
  for( j = 0 ; j < NS ; j++ ) {
    
    const size_t col2 = GSRC.ig[ j ] ;
    
    for( i = 0 ; i < NS ; i++ ) {
      
      const size_t col1 = GSNK.ig[ i ] ;
      
      // cache the color matrices
      fcache = (const __m128d*)fwd.D[j][i].C ;
      bcache = (const __m128d*)bwd.D[col1][col2].C ;
      sum = _mm_setzero_pd( ) ;
      
      // unroll
      #if NC == 3
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[0] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[3] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[6] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[1] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[4] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[7] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[2] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[5] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[8] ) ) ; bcache++ ;
      #elif NC == 2
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[0] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[2] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[1] ) ) ; bcache++ ;
      sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[3] ) ) ; bcache++ ;
      #else 
      size_t c1 , c2 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  sum = _mm_add_pd( sum , SSE2_MUL( *bcache , fcache[ c1 + c2*NC ] ) ) ; bcache++ ;
	}
      }
      #endif

      // switch for the phases (note the implicit minus sign!)
      switch( ( GSNK.g[ i ] + GSRC.g[ col2 ] ) & 3 ) {
      case 0 : gsum = _mm_sub_pd( gsum , sum ) ; break ;
      case 1 : gsum = _mm_add_pd( gsum , _mm_shuffle_pd( sum , SSE_FLIP( sum ) , 1 ) ) ; break ;
      case 2 : gsum = _mm_add_pd( gsum , sum ) ; break ;
      case 3 : gsum = _mm_add_pd( gsum , _mm_shuffle_pd( SSE_FLIP( sum ) , sum , 1 ) ) ; break ;
      }
    }
  }

  double complex s ;
  _mm_store_pd( (void*)&s , gsum ) ;
  return s ;
}

#endif
