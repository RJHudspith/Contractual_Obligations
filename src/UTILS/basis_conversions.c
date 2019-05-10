/**
   @file basis_conversions.c
   @brief conversions between our chiral and the NREL (Dirac-Pauli) basis 
 */
#include "common.h"

#ifdef __SSE2__
  #include <immintrin.h>
#endif

// add four color matrices
static void
add_colors( void *R ,
	    const void *A ,
	    const void *B ,
	    const void *C ,
	    const void *D )
{
#ifdef __SSE2__
  __m128d *pR = (__m128d*)R ;
  const __m128d *pA = (const __m128d*)A , *pB = (const __m128d*)B ;
  const __m128d *pC = (const __m128d*)C , *pD = (const __m128d*)D ;
  register const __m128d onehalf = _mm_setr_pd( 0.5 , 0.5 ) ;
  size_t c ;
  for( c = 0 ; c < NCNC ; c++ ) {
    *pR = _mm_mul_pd( onehalf , _mm_add_pd( _mm_add_pd( *pA , *pB ) ,
					    _mm_add_pd( *pC , *pD ) ) ) ;
    pR++ ; pA++ ; pB++ ; pC++ ; pD++ ;
  }
#else
  double *pR = (double*)R ;
  const double complex *pA = (const double complex*)A ;
  const double complex *pB = (const double complex*)B ;
  const double complex *pC = (const double complex*)C ;
  const double complex *pD = (const double complex*)D ;
  register const double onehalf = 0.5 ;
  size_t c ;
  for( c = 0 ; c < NCNC ; c++ ) {
    *pR = onehalf*( *pA + *pB + *pC + *pD ) ;
    pR++ ; pA++ ; pB++ ; pC++ ; pD++ ;
  }
#endif
}

// do (A+B)-(C+D) for color matrices 
static void
sub_colors( void *R ,
	    const void *A ,
	    const void *B ,
	    const void *C ,
	    const void *D )
{
#ifdef __SSE2__
  __m128d *pR = (__m128d*)R ;
  const __m128d *pA = (const __m128d*)A , *pB = (const __m128d*)B ;
  const __m128d *pC = (const __m128d*)C , *pD = (const __m128d*)D ;
  register const __m128d onehalf = _mm_setr_pd( 0.5 , 0.5 ) ;
  size_t c ;
  for( c = 0 ; c < NCNC ; c++ ) {
    *pR = _mm_mul_pd( onehalf , _mm_sub_pd( _mm_add_pd( *pA , *pB ) ,
					    _mm_add_pd( *pC , *pD ) ) ) ;
    pR++ ; pA++ ; pB++ ; pC++ ; pD++ ;
  }
#else
  double *pR = (double*)res ;
  const double complex *pA = (const double complex*)A ;
  const double complex *pB = (const double complex*)B ;
  const double complex *pC = (const double complex*)C ;
  const double complex *pD = (const double complex*)D ;
  register const double onehalf = 0.5 ;
  size_t c ;
  for( c = 0 ; c < NCNC ; c++ ) {
    *pR = onehalf*( *pA + *pB - ( *pC + *pD ) ) ;
    pR++ ; pA++ ; pB++ ; pC++ ; pD++ ;
  }
#endif
}

// is rotated from the eigenvectors of gamma_t
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  // first row
  add_colors( S->D[0][0].C , P.D[0][0].C , P.D[2][0].C , P.D[0][2].C , P.D[2][2].C ) ;
  add_colors( S->D[0][1].C , P.D[0][1].C , P.D[2][1].C , P.D[0][3].C , P.D[2][3].C ) ;
  sub_colors( S->D[0][2].C , P.D[0][2].C , P.D[2][2].C , P.D[0][0].C , P.D[2][0].C ) ;
  sub_colors( S->D[0][3].C , P.D[0][3].C , P.D[2][3].C , P.D[0][1].C , P.D[2][1].C ) ;
  // second row
  add_colors( S->D[1][0].C , P.D[1][0].C , P.D[3][0].C , P.D[1][2].C , P.D[3][2].C ) ;
  add_colors( S->D[1][1].C , P.D[1][1].C , P.D[3][1].C , P.D[1][3].C , P.D[3][3].C ) ;
  sub_colors( S->D[1][2].C , P.D[1][2].C , P.D[3][2].C , P.D[1][0].C , P.D[3][0].C ) ;
  sub_colors( S->D[1][3].C , P.D[1][3].C , P.D[3][3].C , P.D[1][1].C , P.D[3][1].C ) ;
  // third row
  sub_colors( S->D[2][0].C , P.D[2][0].C , P.D[2][2].C , P.D[0][2].C , P.D[0][0].C ) ;
  sub_colors( S->D[2][1].C , P.D[2][1].C , P.D[2][3].C , P.D[0][3].C , P.D[0][1].C ) ;
  sub_colors( S->D[2][2].C , P.D[0][0].C , P.D[2][2].C , P.D[2][0].C , P.D[0][2].C ) ;
  sub_colors( S->D[2][3].C , P.D[0][1].C , P.D[2][3].C , P.D[2][1].C , P.D[0][3].C ) ;
  // fourth row
  sub_colors( S->D[3][0].C , P.D[3][0].C , P.D[3][2].C , P.D[1][0].C , P.D[1][2].C ) ;
  sub_colors( S->D[3][1].C , P.D[3][1].C , P.D[3][3].C , P.D[1][1].C , P.D[1][3].C ) ;
  sub_colors( S->D[3][2].C , P.D[1][0].C , P.D[3][2].C , P.D[3][0].C , P.D[1][2].C ) ;
  sub_colors( S->D[3][3].C , P.D[1][1].C , P.D[3][3].C , P.D[3][1].C , P.D[1][3].C ) ;	
  return ;
}

// rotate a timeslice
static void
nrel_rotate_slice( struct spinor *S )
{
  size_t site ;
#pragma omp for private(site) 
  for( site = 0 ; site < LCU ; site++ ) {
    chiral_to_nrel( &S[ site ] ) ;
  }
  return ;
}

// rotate if we need to
void 
rotate_offdiag( struct spinor **S ,
		const struct propagator *prop ,
		const size_t Nprops )
{
  // loop all props looking to see if any are non-relativistic
  size_t mu ;
  GLU_bool have_NREL = GLU_FALSE ;
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    switch( prop[mu].basis ) {
    case NREL_FWD  : have_NREL = GLU_TRUE ; break ;
    case NREL_BWD  : have_NREL = GLU_TRUE ; break ;
    case NREL_CORR : have_NREL = GLU_TRUE ; break ;
    case CHIRAL : break ;
    }
  }
  // leave if it is all chiral
  if( have_NREL == GLU_FALSE ) return ;

  // loop back through the list rotating any chiral props
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    switch( prop[mu].basis ) {
    case NREL_FWD  : break ;
    case NREL_BWD  : break ;
    case NREL_CORR : break ;
    case CHIRAL :
      nrel_rotate_slice( S[ mu ] ) ;
      break ;
    }
  }
  
  return ;
}
