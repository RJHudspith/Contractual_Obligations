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
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      dagger_gauge( (double complex*)adj -> D[d2][d1].C ,
		    (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// returns spin-color trace : Tr[ A B ]
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B )
{
  int d1 , d2 ;
  register double complex sum = 0.0 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      sum += colortrace_prod( (const double complex*)A.D[d1][d2].C ,
			      (const double complex*)B.D[d2][d1].C ) ;
    }
  }
  return sum ;
}

// computes G5 ( adj( S ) ) G5
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
  int i , j ;
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col = GAMMA.ig[i] ;
    double complex fac = 1.0 ;
    // switch for the phases
    switch( GAMMA.g[i] ) {
    case 0 : fac =  1 ; break ;
    case 1 : fac =  I ; break ;
    case 2 : fac = -1 ; break ;
    case 3 : fac = -I ; break ;
    }
    // multiply out
    for( j = 0 ; j < NS ; j++ ) {
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)res -> D[col][j].C ) ;
    } 
    //
  }
  //
  *res = tmp ;
  return ;
}

// multiply a spinor on the right with a gamma matrix
void
gamma_mul_r( struct spinor *__restrict res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  int i , j ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const int col = GAMMA.ig[j] ;
    double complex fac = 1.0 ;
    switch( GAMMA.g[ col ] ) {
    case 0 : fac =  1 ; break ;
    case 1 : fac =  I ; break ;
    case 2 : fac = -1 ; break ;
    case 3 : fac = -I ; break ;
    }
    // and copy it in
    for( i = 0 ; i < NS ; i++ ) {
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)res -> D[i][col].C ) ;
    } 
    //
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
  double complex fac = 1.0 ;
  int i , j ;
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col1 = GLEFT.ig[ i ] ;
    for( j = 0 ; j < NS ; j++ ) {
      const int col2 = GRIGHT.ig[ j ] ;
      // switch for the phases
      switch( ( GLEFT.g[ i ] + GRIGHT.g[ col2 ] ) & 3 ) {
      case 0 : fac =  1 ; break ;
      case 1 : fac =  I ; break ;
      case 2 : fac = -1 ; break ;
      case 3 : fac = -I ; break ;
      }
      // multiply out
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)S -> D[col1][col2].C ) ;
    }
  }
  //
  *S = tmp ;
  return ;
}

static inline __m128d
SSE2_MULCONJ( const __m128d a , const __m128d b )
{
  return _mm_add_pd( _mm_mul_pd( _mm_unpacklo_pd( a , a ) , b ) ,                               // re(a)*re(b) , re(a)*im(b)
		     _mm_mul_pd( _mm_unpackhi_pd( a , -a ) , _mm_shuffle_pd( b , b , 1 ) ) ) ;  // im(a)*im(b) , -re(b)*im(a)
}

// meson contraction code computes Tr[ GSNK ( G5 bwd G5 )^{\dagger} GSRC ( fwd ) ]
double complex
meson_contract( const struct gamma GSNK ,		
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 )
{
  register __m128d gsum = _mm_setzero_pd( ) ;
  const __m128d *d2 = (const __m128d*)fwd.D ;

  int i , j , col2 , col1 ;
  for( j = 0 ; j < NS ; j++ ) {
    
    col2 = GSRC.ig[ G5.ig[ j ] ] ; 

    // loop columns
    for( i = 0 ; i < NS ; i++ ) {

      col1 = G5.ig[ GSNK.ig[ i ] ] ;
      
      // sums in double to avoid complex multiply
      const __m128d *d1 = (const __m128d*)bwd.D[col2][col1].C ;

      #if NC == 3
      register __m128d sum = SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      #else 
      register __m128d sum = _mm_setzero_pd( ) ;
      int c1c2 ;
      for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
	sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      }
      #endif

      // switch for the phases
      switch( ( GSNK.g[ i ] + G5.g[ col1 ] + G5.g[ col2 ] + GSRC.g[ col2 ] ) & 3 ) {
      case 0 : gsum += sum ; break ;
      case 1 : gsum += _mm_shuffle_pd( -sum , sum , 1 ) ; break ;
      case 2 : gsum += -sum ; break ;
      case 3 : gsum += _mm_shuffle_pd( sum , -sum , 1 ) ; break ;
      }

      // and we are done
    }
  }

  // cast through the void
  double complex s ;
  _mm_store_pd( (void*)&s , gsum ) ;
  return s ;
}

// meson contraction code computes Tr[ GSNK ( bwd ) GSRC ( fwd ) ]
double complex
simple_meson_contract( const struct gamma GSNK ,		
		       const struct spinor bwd , 
		       const struct gamma GSRC ,
		       const struct spinor fwd )
{
  register __m128d gsum = _mm_setzero_pd( ) ;
  const __m128d *d2 = (const __m128d*)fwd.D ;

  int i , j ;
  for( i = 0 ; i < NS ; i++ ) {

    const int col1 = GSNK.ig[ i ] ;

    // loop columns
    for( j = 0 ; j < NS ; j++ ) {

      const int col2 = GSRC.ig[ j ] ;

      // sums in double to avoid complex multiply
      const __m128d *d1 = (const __m128d*)bwd.D[col2][col1].C ;

      #if NC == 3
      register __m128d sum = SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      #else 
      register __m128d sum = _mm_setzero_pd( ) ;
      int c1c2 ;
      for( c1c2 = 0 ; c1c2 < NCNC ; c1c2++ ) {
	sum += SSE2_MULCONJ( *d2 , *d1 ) ; d1 ++ ; d2 ++ ;
      }
      #endif

      // switch for the phases
      switch( ( GSNK.g[ i ] + GSRC.g[ col2 ] ) & 3 ) {
      case 0 : gsum += sum ; break ;
      case 1 : gsum += _mm_shuffle_pd( -sum , sum , 1 ) ; break ;
      case 2 : gsum += -sum ; break ;
      case 3 : gsum += _mm_shuffle_pd( sum , -sum , 1 ) ; break ;
      }
      //
    }
  }

  // cast through the void
  double complex s ;
  _mm_store_pd( (void*)&s , gsum ) ;
  return s ;
}

#endif
