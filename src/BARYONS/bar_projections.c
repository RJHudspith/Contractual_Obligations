/**
   @file bar_projections.c
   @brief spin-1/2 and spin-3/2 projections
 */
#include "common.h"

#include "gammas.h"             // gamma_mmul()
#include "spinmatrix_ops.h"     // matrix operations

// compute psq
void
compute_p_psq( double p[ ND ] ,
	       double *p2 ,
	       const struct veclist momentum )
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    p[ mu ] = sin( momentum.MOM[ mu ] * Latt.twiddles[ mu ] ) ;
    *p2 += p[ mu ] * p[ mu ] ;
  }
  return ;
}

// is the factor proj = ( pslash \gamma_i p_j )
static void
momentum_factor_left( double complex *proj ,
		      const double complex *pslash ,
		      const struct gamma *GAMMA ,
		      const double p[ ND ] ,
		      const size_t i ,
		      const size_t j )
{
  // temporary space
  double complex *tmp1 = NULL , *tmp2 = NULL ;
  if( corr_malloc( (void**)&tmp1 , 16 , NSNS*sizeof( double complex ) ) != 0 ||
      corr_malloc( (void**)&tmp2 , 16 , NSNS*sizeof( double complex ) ) != 0 ) {
    goto memfree ;
  }
  // set proj to zero
  zero_spinmatrix( proj ) ;
  // no use in adding zeros
  if( fabs( p[ j ] ) > 0 ) {
    // set tmp1 to be diagonal with p[ j ]
    zero_spinmatrix( tmp1 ) ;
    size_t d1  ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      tmp1[ d1 * ( NS + 1 ) ] = p[ j ] ;
    }
    // left multiply by gamma_i
    spinmatrix_gamma( tmp1 , GAMMA[ i ] ) ;
    // left multiply by pslash
    spinmatrix_multiply( tmp2 , pslash , tmp1 ) ;
    // add into D
    atomic_add_spinmatrices( proj , tmp2 ) ;
  }
 memfree :
  free( tmp1 ) ; free( tmp2 ) ;
  return ;
}

// computes the right hand factor proj = ( p_i \gamma_j pslash )
static void
momentum_factor_right( double complex *proj ,
		       const double complex *pslash ,
		       const struct gamma *GAMMA ,
		       const double p[ ND ] ,
		       const size_t i ,
		       const size_t j )
{
  // temporary space
  double complex *tmp1 = NULL , *tmp2 = NULL ;
  if( corr_malloc( (void**)&tmp1 , 16 , NSNS*sizeof( double complex ) ) != 0 ||
      corr_malloc( (void**)&tmp2 , 16 , NSNS*sizeof( double complex ) ) != 0 ) {
    goto memfree ;
  }
  // set proj to zero
  zero_spinmatrix( proj ) ;
  // no use in adding zeros
  if( fabs( p[ i ] ) > 0 ) {
    // set tmp1 to be diagonal with p[ i ]
    zero_spinmatrix( tmp1 ) ;
    size_t d1 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      tmp1[ d1 * ( NS + 1 ) ] = p[ i ] ;
    }
    // right multiply by gamma_j
    spinmatrix_gamma( tmp1 , GAMMA[ j ] ) ;
    // right multiply by pslash
    spinmatrix_multiply( tmp2 , tmp1 , pslash ) ;
    // add into D
    atomic_add_spinmatrices( proj , tmp2 ) ;
  }
 memfree :
  free( tmp1 ) ; free( tmp2 ) ;
  return ;
}

// p00 projection matrix
void
P00( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  zero_spinmatrix( proj ) ;
  return ;
}

// p11 projection matrix
void
P11( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  // set projection to the identity
  identity_spinmatrix( proj ) ;

  // zero momentum bit first
  struct gamma Gij ; 
  gamma_mmul( &Gij , GAMMA[ i ] , GAMMA[ j ] ) ;
  spinmatrix_gamma( proj , Gij ) ;
  spinmatrix_mulconst( proj , 1.0 / ( ND - 1 ) ) ;

  // compute p2
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;

  // add the momentum factors
  if( p2 > 0 ) {
    double complex *left = NULL , *right = NULL , *pslash = NULL ;
    size_t d ;
    if( corr_malloc( (void**)&left , 16 , NSNS*sizeof( double complex ) ) != 0 ||
	corr_malloc( (void**)&right , 16 , NSNS*sizeof( double complex ) ) != 0 || 
	corr_malloc( (void**)&pslash , 16 , NSNS*sizeof( double complex ) ) != 0  ) {
      goto memfree ;
    }
    // subtract on the diagonal
    for( d = 0 ; d < NS ; d++ ) {
      proj[ d*( NS+1 ) ] -= ( p[ i ] * p[ j ] ) / p2 ;
    }
    // add the left hand side and right hand side
    compute_pslash( pslash , GAMMA , p ) ;
    // compute left and right hand side
    momentum_factor_left( left , pslash , GAMMA , p , i , j ) ;
    momentum_factor_right( right , pslash , GAMMA , p , i , j ) ;
    atomic_add_spinmatrices( left , right ) ;
    spinmatrix_mulconst( left , 1.0 / ( p2 * ( ND - 1 ) ) ) ;
    // and add it back in
    atomic_add_spinmatrices( proj , left ) ;
  memfree:
    free( left ) ; free( right ) ; free( pslash ) ;
  }
  return ;
}

// p12 projection matrix
void
P12( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  // set projection to zero
  zero_spinmatrix( proj ) ;

  // compute p2
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;

  if( p2 > 0 ) {
    double complex *pslash = NULL ;
    if( corr_malloc( (void**)&pslash , 16 , NSNS*sizeof( double complex ) ) != 0 ) {
      goto memfree ;
    }
    compute_pslash( pslash , GAMMA , p ) ;
    // compute pslash.\gamma_i p_j
    momentum_factor_left( proj , pslash , GAMMA , p , i , j ) ;
    // subtract on the diagonal
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      proj[ d*( NS+1 ) ] -= p[ i ] * p[ j ] ;
    }
    // multiply by -1 / ( sqrt(3) * p^2 )
    spinmatrix_mulconst( proj , -1.0 / ( sqrt( 3 ) * p2 ) ) ;
  memfree : 
    free( pslash ) ;
  }
  return ;
}

// p12 projection matrix
void
P21( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  // set to zero
  zero_spinmatrix( proj ) ;

  // compute p2
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;

  if( p2 > 0 ) {
    // compute pslash
    double complex *pslash = NULL ;
    if( corr_malloc( (void**)&pslash , 16 , NSNS*sizeof( double complex ) ) != 0 ) {
      goto memfree ;
    }
    compute_pslash( pslash , GAMMA , p ) ;
    // compute pslash.\gamma_j.p_i
    momentum_factor_left( proj , pslash , GAMMA , p , j , i ) ;
    // subtract the on the diagonal
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      proj[ d*( NS+1 ) ] -= p[ i ] * p[ j ] ;
    }
    // multiply by -1 / ( sqrt(3) * p^2 )
    spinmatrix_mulconst( proj , 1.0 / ( sqrt( 3 ) * p2 ) ) ;
  memfree :
    free( pslash ) ;
  }
  return ;
}

// p22 projection matrix
void
P22( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  // set to zero
  zero_spinmatrix( proj ) ;

  // compute p2
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;

  if( p2 > 0 && fabs( p[ i ] ) > 0 && fabs( p[ j ] ) > 0 ) {
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      proj[ d * ( NS + 1 ) ] = p[ i ] * p[ j ] / p2 ;
    }
  }
  return ;
}

// p32 projection matrix
void
P32( double complex *proj ,
     const size_t i ,
     const size_t j ,
     const struct gamma *GAMMA ,
     const struct veclist *momentum ,
     const size_t pidx )
{
  // set proj to the identity
  identity_spinmatrix( proj ) ;

  // zero momentum bit first
  struct gamma Gij ; 
  gamma_mmul( &Gij , GAMMA[ i ] , GAMMA[ j ] ) ;
  spinmatrix_gamma( proj , Gij ) ;
  spinmatrix_mulconst( proj , -1.0 / ( ND - 1 ) ) ;

  // compute p2
  double p2 = 0.0 , p[ ND ] ;
  compute_p_psq( p , &p2 , momentum[ pidx ] ) ;

  // add the momentum factors
  if( p2 > 0 ) {
    // add the left hand side and right hand side
    double complex *left = NULL , *right = NULL , *pslash = NULL ;
    if( corr_malloc( (void**)&left , 16 , NSNS*sizeof( double complex ) ) != 0 ||
	corr_malloc( (void**)&right , 16 , NSNS*sizeof( double complex ) ) != 0 || 
	corr_malloc( (void**)&pslash , 16 , NSNS*sizeof( double complex ) ) != 0  ) {
      goto memfree ;
    }
    compute_pslash( pslash , GAMMA , p ) ;
    // compute left and right hand side
    momentum_factor_left( left , pslash , GAMMA , p , i , j ) ;
    momentum_factor_right( right , pslash , GAMMA , p , i , j ) ;
    atomic_add_spinmatrices( left , right ) ;
    spinmatrix_mulconst( left , -1.0 / ( p2 * ( ND - 1 ) ) ) ;
    // and add it back in
    atomic_add_spinmatrices( proj , left ) ;
  memfree : 
    free( left ) ; free( right ) ; free( pslash ) ; 
  }

  // delta function g_ij
  if( i == j ) {
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      proj[ d * ( NS + 1 ) ] += 1.0 ;
    }
  }
  return ;
}

// select matrix G_{ik} with dirac indices i (row) k (cols)
static void
set_Gik( double complex *Gik , 
	 const struct mcorr **corr ,
	 const size_t i ,
	 const size_t k ,
	 const size_t pidx ,
	 const size_t t )
{
  size_t d ;
#if NS == 4
  for( d = 0 ; d < NS ; d++ ) {
    Gik[ 0 + d * NS ] = corr[ k + i*B_CHANNELS ][ 0 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 1 + d * NS ] = corr[ k + i*B_CHANNELS ][ 1 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 2 + d * NS ] = corr[ k + i*B_CHANNELS ][ 2 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 3 + d * NS ] = corr[ k + i*B_CHANNELS ][ 3 + d * NS ].mom[ pidx ].C[ t ] ;
  }
#else
  for( d = 0 ; d < NSNS ; d++ ) {
    Gik[ d ] = corr[ k + i*B_CHANNELS ][ d ].mom[ pidx ].C[ t ] ;
  }
#endif
}

// project out a spinmatrix into D
// performs G_{ik} = P_{ij} G_{jk}
int
spinproj( double complex *Gik ,
	  void (*p)( double complex *proj ,
		     const size_t i ,
		     const size_t j ,
		     const struct gamma *GAMMA ,
		     const struct veclist *momentum ,
		     const size_t pidx ) ,
	  const struct mcorr **corr ,
	  const size_t i ,
	  const size_t k ,
	  const size_t t ,
	  const struct gamma *GAMMA ,
	  const struct veclist *momentum ,
	  const size_t pidx )
{
  // set sum to zero
  zero_spinmatrix( Gik ) ;

  // absolute value of momenta
  size_t mu , sum = 0 , j ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    sum += abs( momentum[ pidx ].MOM[ mu ] ) ;
  }

  // flag for success or failure of the routine
  size_t flag = SUCCESS ; 

  // temporary storage for the projector
  double complex *pij = NULL , *tmp = NULL , *Gjk = NULL ;
  if( corr_malloc( (void**)&pij , 16 , NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }
  if( corr_malloc( (void**)&tmp , 16 , NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }
  if( corr_malloc( (void**)&Gjk , 16 , NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }

  // perform contraction
  for( j = 0 ; j < ND ; j++ ) {
    set_Gik( Gjk , corr , j , k , pidx , t ) ;
    p( pij , i , j , GAMMA , momentum , pidx ) ;
    spinmatrix_multiply( tmp , pij , Gjk ) ;
    atomic_add_spinmatrices( Gik , tmp ) ;
  }

 memfree :
    
  free( tmp ) ;
  free( pij ) ;
  free( Gjk ) ;

  return flag ;
}

// project out a specific spin
static void
spin_project( double complex *Gik , 
	      const struct mcorr **corr ,
	      const struct gamma *GAMMA ,
	      const struct veclist *momentum ,
	      const size_t i ,
	      const size_t k ,
	      const size_t p ,
	      const size_t t ,
	      const spinhalf projection ) 
{
  // spin projections are for spatial indices only!!
  if( i > ( ND-2 ) || k > ( ND-2 ) ) {
    set_Gik( Gik , corr , i , k , p , t ) ;
    return ;
  }
  // do the spin projections
  switch( projection ) {
  case OneHalf_11 :
    spinproj( Gik , P11 , corr , i , k , t , GAMMA , momentum , p ) ;
    break ;
  case OneHalf_12 :
    spinproj( Gik , P12 , corr , i , k , t , GAMMA , momentum , p ) ;
    break ;
  case OneHalf_21 :
    spinproj( Gik , P21 , corr , i , k , t , GAMMA , momentum , p ) ;
    break ;
  case OneHalf_22 :
    spinproj( Gik , P22 , corr , i , k , t , GAMMA , momentum , p ) ;
    break ;
  case ThreeHalf :
    spinproj( Gik , P32 , corr , i , k , t , GAMMA , momentum , p ) ;
    break ;
  case NONE :
    set_Gik( Gik , corr , i , k , p , t ) ;
    break ;
  }
  return ;
}

// returns the correlator at some momentum
double complex*
baryon_project( const struct mcorr **corr ,
		const struct gamma *GAMMA ,
		const struct veclist *momentum ,
		const size_t GSRC ,
		const size_t GSNK ,
		const size_t p ,
		const bprojection parity_proj , 
		const spinhalf spin_proj ) 
{
  // D is a temporary spinmatrix
  double complex *D = NULL , *result = NULL ;

  // allocate temporaries
  corr_malloc( (void**)&D , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&result , 16 , LT * sizeof( double complex ) ) ;

  // build these just in case
  struct gamma g0g5 ; gamma_mmul( &g0g5 , GAMMA[ 0 ] , GAMMA[ 5 ] ) ;
  struct gamma g3g0g5 ; gamma_mmul( &g3g0g5 , GAMMA[ 3 ] , g0g5 ) ;

  // loop times 
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {

    // set open dirac matrix doing the spin projection if desired
    spin_project( D , corr , GAMMA , momentum , 
		  GSRC , GSNK , p , t , spin_proj ) ; 

    switch( parity_proj ) {
    case L0 : // does 1/4( g_I + g_3 - i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) + 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) - 
			     I * gammaspinmatrix_trace( g0g5 , D ) - 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L1 : // does 1/4( g_I - g_3 + i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) +
			     I * gammaspinmatrix_trace( g0g5 , D ) - 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L2 : // does 1/4( g_I - g_3 - i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) -
			     I * gammaspinmatrix_trace( g0g5 , D ) + 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L3 : // does 1/4( g_I + g_3 + i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) + 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) +
			     I * gammaspinmatrix_trace( g0g5 , D ) + 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L4 : // does 1/2 ( g_I + g_3 )
      result[ t ] = 0.5 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) +
			    gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) ) ;
      break ;
    case L5 : // does 1/2 ( g_I - g_3 )
      result[ t ] = 0.5 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			    gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) ) ;
      break ;
    }
  }

  free( D ) ;
  return result ;
}
