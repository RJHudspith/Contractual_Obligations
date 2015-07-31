/**
   @file baryon_projections.c
   @brief projections available for our baryons

   Taken from the papers 
 */
#include "common.h"

#include "baryon_projections.h" // for the enums
#include "gammas.h"             // gamma_mmul()
#include "spinmatrix_ops.h"     // matrix operations

// is the factor proj = ( pslash \gamma_i p_j )
static void
momentum_factor_left( double complex *proj ,
		      const double complex *pslash ,
		      const struct gamma *GAMMA ,
		      const size_t DIMS ,
		      const double p[ DIMS ] ,
		      const size_t i ,
		      const size_t j )
{
  double complex tmp2[ NSNS ] , tmp1[ NSNS ] ;
  size_t d1  ;
  // set proj to zero
  zero_spinmatrix( proj ) ;
  // no use in adding zeros
  if( p[ j ] > 0 ) {
    // set tmp1 to be diagonal with p[ j ]
    zero_spinmatrix( tmp1 ) ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      tmp1[ d1 * ( NS + 1 ) ] = p[ j ] ;
    }
    // left multiply by gamma_i
    gamma_spinmatrix( tmp1 , GAMMA[ i ] ) ;
    // left multiply by pslash
    spinmatrix_multiply( tmp2 , pslash , tmp1 ) ;
    // add into D
    atomic_add_spinmatrices( proj , tmp2 ) ;
  }
  return ;
}

// computes the right hand factor proj = ( p_i \gamma_j pslash )
static void
momentum_factor_right( double complex *proj ,
		       const double complex *pslash ,
		       const struct gamma *GAMMA ,
		       const size_t DIMS ,
		       const double p[ DIMS ] ,
		       const size_t i ,
		       const size_t j )
{
  double complex tmp2[ NSNS ] , tmp1[ NSNS ] ;
  size_t d1 ;
  // set proj to zero
  zero_spinmatrix( proj ) ;
  // no use in adding zeros
  if( p[ i ] > 0 ) {
    // set tmp1 to be diagonal with p[ i ]
    zero_spinmatrix( tmp1 ) ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      tmp1[ d1 * ( NS + 1 ) ] = p[ j ] ;
    }
    // right multiply by gamma_j
    spinmatrix_gamma( tmp1 , GAMMA[ j ] ) ;
    // right multiply by pslash
    spinmatrix_multiply( tmp2 , tmp1 , pslash ) ;
    // add into D
    atomic_add_spinmatrices( proj , tmp2 ) ;
  }
  return ;
}

// select matrix with dirac indices j (row) k (cols)
static void
set_Ojk( double complex *Ojk , 
	 const struct mcorr **corr ,
	 const size_t j ,
	 const size_t k ,
	 const size_t pidx ,
	 const size_t t ,
	 const double factor )
{
  size_t d ;
#if NS == 4
  for( d = 0 ; d < NS ; d++ ) {
    Ojk[ 0 + d * NS ] = corr[ k + j*B_CHANNELS ][ 0 + d * NS ].mom[ pidx ].C[ t ] * factor ;
    Ojk[ 1 + d * NS ] = corr[ k + j*B_CHANNELS ][ 1 + d * NS ].mom[ pidx ].C[ t ] * factor ;
    Ojk[ 2 + d * NS ] = corr[ k + j*B_CHANNELS ][ 2 + d * NS ].mom[ pidx ].C[ t ] * factor ;
    Ojk[ 3 + d * NS ] = corr[ k + j*B_CHANNELS ][ 3 + d * NS ].mom[ pidx ].C[ t ] * factor ;
  }
#else
  for( d = 0 ; d < NSNS ; d++ ) {
    Ojk[ d ] = corr[ k + j*B_CHANNELS ][ d ].mom[ pidx ].C[ t ] * factor ;
  }
#endif
}

// compute the momenta
static void
compute_p_psq( const size_t DIMS ,
	       double p[ DIMS ] ,
	       double *p2 ,
	       const struct veclist momentum )
{
  size_t mu ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    p[ mu ] = sin( momentum.MOM[mu] * Latt.twiddles[mu] ) ;
    *p2 += p[ mu ] * p[ mu ] ;
  }
  return ;
}

// spin 1/2 projection (P_11)_{ij} O_{jk} = ( \gamma_i \gamma_j - p_i p_j / ( 3 * p^2 ) + 
//                                            1/(3p^2)( pslash \gamma_i p_j + p_i \gamma_j pslash ) ) O_{jk} summed into D
static void
spinhalf11_project( double complex *D , 
		    const struct mcorr **corr ,
		    const struct gamma *GAMMA ,
		    const struct veclist *momentum ,
		    const size_t i ,
		    const size_t k ,
		    const size_t pidx ,
		    const size_t t )
{
  double complex *tmp1 = NULL , *pslash = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&pslash , 16 , NSNS * sizeof( double complex ) ) ;
  // set result to 0
  zero_spinmatrix( D ) ;
  // compute p2
  double p2 = 0.0 , p[ ND - 1 ] ;
  compute_p_psq( ND-1 , p , &p2 , momentum[ pidx ] ) ;
  // compute pslash
  compute_pslash( pslash , GAMMA , ND-1 , p ) ;
  size_t j ;
  // perform 1/3 \gamma_i \gamma_j O_j
  for( j = 0 ; j < ND-1 ; j++ ) {
    // set it up
    set_Ojk( tmp1 , corr , j , k , pidx , t , 1.0 / 3.0 ) ;
    struct gamma gij ; gamma_mmul( &gij , GAMMA[i] , GAMMA[j] ) ;
    gamma_spinmatrix( tmp1 , gij ) ;
    atomic_add_spinmatrices( D , tmp1 ) ;
    // D += ( -p_i p_j / p^2 O_j )
    if( p2 > 0 ) {
      if( p[ i ] > 0 && p[ j ] > 0 ) {
	// set tmp1 to be O_{jk}
	set_Ojk( tmp1 , corr , j , k , pidx , t , 1.0 / 3.0 ) ;
	spinmatrix_mulconst( tmp1 , -p[ i ] * p[ j ] / p2 ) ;
	atomic_add_spinmatrices( D , tmp1 ) ;
      }
      // set tmp1 to be +1/(3p^2) O_j
      set_Ojk( tmp1 , corr , j , k , pidx , t , 1.0 / ( 3.0 * p2 ) ) ;
      double complex proj[ NSNS ] , res[ NSNS ] ;
      // computes the projector and puts it in "proj"
      momentum_factor_left( proj , pslash , GAMMA , ND-1 , p , i , j ) ;
      momentum_factor_right( res , pslash , GAMMA , ND-1 , p , i , j ) ;
      atomic_add_spinmatrices( proj , res ) ;
      // multiply res = proj_{ij} * O_j
      spinmatrix_multiply( res , proj , tmp1 ) ;
      atomic_add_spinmatrices( D , res ) ;
    }
  }
  free( pslash ) ;
  free( tmp1 ) ;
  return ;
}

// spin 1/2 projection (P_12)_{ij} O_{jk} = ( 1/( \sqrt(3)p^2 ) ( p_i p_j - pslash \gamma_i p_j ) O_{jk} summed into D
static void
spinhalf12_project( double complex *D , 
		    const struct mcorr **corr ,
		    const struct gamma *GAMMA ,
		    const struct veclist *momentum ,
		    const size_t i ,
		    const size_t k ,
		    const size_t pidx ,
		    const size_t t )
{
  double complex *tmp1 = NULL , *pslash = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&pslash , 16 , NSNS * sizeof( double complex ) ) ;
  // set result to 0
  zero_spinmatrix( D ) ;
  // compute p2
  double p2 = 0.0 , p[ ND - 1 ] ;
  compute_p_psq( ND-1 , p , &p2 , momentum[ pidx ] ) ;
  // compute pslash
  compute_pslash( pslash , GAMMA , ND-1 , p ) ;
  // compute the projector
  if( p2 > 0 ) {
    size_t j ;
    for( j = 0 ; j < ND-1 ; j++ ) {
      // D += ( p_i p_j ) O_j
      if( p[ i ] > 0 && p[ j ] > 0 ) {
	set_Ojk( tmp1 , corr , j , k , pidx , t , 1.0 / ( sqrt( 3.0 ) * p2 ) ) ;
	spinmatrix_mulconst( tmp1 , p[ i ] * p[ j ] / p2 ) ;
	atomic_add_spinmatrices( D , tmp1 ) ;
      }
      // D -= ( pslash gamma_i p_j ) O_j
      if( p[ j ] > 0 ) {
	// compute O_{jk} * ( 1 / ( sqrt(3) * p^2 ) )
	set_Ojk( tmp1 , corr , j , k , pidx , t , -1.0 / ( sqrt( 3.0 ) * p2 ) ) ;
	double complex proj[ NSNS ] , res[ NSNS ] ;
	// computes the momentum projection and puts in "proj"
	momentum_factor_left( proj , pslash , GAMMA , ND-1 , p , i , j ) ;
	spinmatrix_multiply( res , proj , tmp1 ) ;
	atomic_add_spinmatrices( D , res ) ;
      }
      // and we are done
    }
  }
  free( pslash ) ;
  free( tmp1 ) ;
  return ;
}

// spin 1/2 projection (P_12)_{ij} O_{jk} = ( 1/( \sqrt(3)p^2 ) ( p_i p_j - pslash \gamma_i p_j ) O_j summed into D
static void
spinhalf21_project( double complex *D , 
		    const struct mcorr **corr ,
		    const struct gamma *GAMMA ,
		    const struct veclist *momentum ,
		    const size_t i ,
		    const size_t k ,
		    const size_t pidx ,
		    const size_t t )
{
  double complex *tmp1 = NULL , *pslash = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&pslash , 16 , NSNS * sizeof( double complex ) ) ;
  // set result to 0
  zero_spinmatrix( D ) ;
  // compute p2
  double p2 = 0.0 , p[ ND - 1 ] ;
  compute_p_psq( ND-1 , p , &p2 , momentum[ pidx ] ) ;
  // compute pslash
  compute_pslash( pslash , GAMMA , ND-1 , p ) ;
  // compute the projector
  if( p2 > 0 ) {
    size_t j ;
    for( j = 0 ; j < ND-1 ; j++ ) {
      // D -= ( p_i p_j ) O_j
      if( p[ i ] > 0 && p[ j ] > 0 ) {
	set_Ojk( tmp1 , corr , j , k , pidx , t , -1.0 / ( sqrt( 3.0 ) * p2 ) ) ;
	spinmatrix_mulconst( tmp1 , p[ i ] * p[ j ] / p2 ) ;
	atomic_add_spinmatrices( D , tmp1 ) ;
      }
      // D += ( p_i gamma_j pslash ) O_j
      if( p[ j ] > 0 ) {
	// compute O_j * ( 1 / ( sqrt(3) * p^2 ) )
	set_Ojk( tmp1 , corr , j , k , pidx , t , 1.0 / ( sqrt( 3.0 ) * p2 ) ) ;
	double complex proj[ NSNS ] , res[ NSNS ] ;
	// computes the momentum projection and puts in "proj"
	momentum_factor_right( proj , pslash , GAMMA , ND-1 , p , i , j ) ;
	spinmatrix_multiply( res , proj , tmp1 ) ;
	atomic_add_spinmatrices( D , res ) ;
      }
      // and we are done
    }
  }
  free( pslash ) ;
  free( tmp1 ) ;
  return ;
}

// spin 1/2 projection (P_12)_{ij} G_{jk} = ( p_i p_j / ( p^2 ) ) O_j summed into D
static void
spinhalf22_project( double complex *D , 
		    const struct mcorr **corr ,
		    const struct gamma *GAMMA ,
		    const struct veclist *momentum ,
		    const size_t i ,
		    const size_t k ,
		    const size_t pidx ,
		    const size_t t )
{
  double complex *tmp1 = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  // set result to 0
  zero_spinmatrix( D ) ;
  // compute p2
  double p2 = 0.0 , p[ ND - 1 ] ;
  compute_p_psq( ND-1 , p , &p2 , momentum[ pidx ] ) ;
  // perform ( p_i p_j ) / ( p^2 )
  if( p2 > 0 ) {
    size_t j ;
    for( j = 0 ; j < ND-1 ; j++ ) {
      set_Ojk( tmp1 , corr , j , k , pidx , t , p[ i ] * p[ j ] / p2 ) ;
      atomic_add_spinmatrices( D , tmp1 ) ;
    }
  }
  free( tmp1 ) ;
  return ;
}

// spin 3/2 projection? summed into D
// is ( delta_{ij} - \gamma_i \gamma_j / 3 - 
//      1/(3p^2) ( pslash \gamma_i p_j + p_i \gamma_j pslash ) )
static void
spinthreehalf_project( double complex *D , 
		       const struct mcorr **corr ,
		       const struct gamma *GAMMA ,
		       const struct veclist *momentum ,
		       const size_t i ,
		       const size_t k ,
		       const size_t pidx ,
		       const size_t t )
{
  double complex *tmp1 = NULL , *pslash = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&pslash , 16 , NSNS * sizeof( double complex ) ) ;
  // set result to 0
  zero_spinmatrix( D ) ;
  // compute p2
  double p2 = 0.0 , p[ ND - 1 ] ;
  compute_p_psq( ND-1 , p , &p2 , momentum[ pidx ] ) ;
  // compute pslash
  compute_pslash( pslash , GAMMA , ND-1 , p ) ;
  size_t j ;
  // perform delta_ij - 1/3 \gamma_i \gamma_j G_{ji}
  for( j = 0 ; j < ND-1 ; j++ ) {
    set_Ojk( tmp1 , corr , j , k , pidx , t , -1.0 / 3.0 ) ;
    if( i == j ) atomic_add_spinmatrices( D , tmp1 ) ;
    struct gamma gij ; gamma_mmul( &gij , GAMMA[i] , GAMMA[j] ) ;
    gamma_spinmatrix( tmp1 , gij ) ;
    atomic_add_spinmatrices( D , tmp1 ) ;
    // if p2 > 0 we do a momentum projection
    if( p2 > 0 ) {
      // set tmp1 to be -1/(3p^2) O_i
      set_Ojk( tmp1 , corr , j , k , pidx , t , -1.0 / ( 3.0 * p2 ) ) ;
      double complex proj[ NSNS ] , res[ NSNS ] ;
      // computes the momentum projection and puts in "proj"
      momentum_factor_left( proj , pslash , GAMMA , ND-1 , p , i , j ) ;
      momentum_factor_right( res , pslash , GAMMA , ND-1 , p , i , j ) ;
      atomic_add_spinmatrices( proj , res ) ;
      // and multiply
      spinmatrix_multiply( res , proj , tmp1 ) ;
      atomic_add_spinmatrices( D , res ) ;
    }
  }
  free( pslash ) ;
  free( tmp1 ) ;
  return ;
}

// projection is P_{ik} G_{kl} = G_{il}
void
spin_project( double complex *Oi , 
	      const struct mcorr **corr ,
	      const struct gamma *GAMMA ,
	      const struct veclist *momentum ,
	      const size_t GSRC ,
	      const size_t GSNK ,
	      const size_t p ,
	      const size_t t ,
	      const spinhalf projection ) 
{
  if( GSRC > 2 || GSNK > 2 ) {
    size_t d ;
    for( d = 0 ; d < NSNS ; d++ ) {
      Oi[ d ] = corr[ GSNK + B_CHANNELS * GSRC ][ d ].mom[ p ].C[ t ] ;
    }
    return ;
  }
  switch( projection ) {
  case OneHalf_11 :
    spinhalf11_project( Oi , corr , GAMMA , momentum ,
			GSRC , GSNK , p , t ) ;
    break ;
  case OneHalf_12 :
    spinhalf12_project( Oi , corr , GAMMA , momentum ,
			GSRC , GSNK , p , t ) ;
    break ;
  case OneHalf_21 :
    spinhalf21_project( Oi , corr , GAMMA , momentum ,
			GSRC , GSNK , p , t ) ;
    break ;
  case OneHalf_22 :
    spinhalf22_project( Oi , corr , GAMMA , momentum ,
			GSRC , GSNK , p , t ) ;
    break ;
  case ThreeHalf :
    spinthreehalf_project( Oi , corr , GAMMA , momentum , 
			   GSRC , GSNK , p , t ) ;
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
		const bprojection projection ) 
{
  // D is a temporary spinmatrix
  double complex *D = NULL , *result = NULL ;

  corr_malloc( (void**)&D , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&result , 16 , LT * sizeof( double complex ) ) ;

  // build these just in case
  struct gamma g0g5 ; gamma_mmul( &g0g5 , GAMMA[ 0 ] , GAMMA[ 5 ] ) ;
  struct gamma g3g0g5 ; gamma_mmul( &g3g0g5 , GAMMA[ 3 ] , g0g5 ) ;

  // loop times
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {
    
    // set D
    set_Ojk( D , corr , GSRC , GSNK , p , t , 1.0 ) ;

    switch( projection ) {
    case L0 : // does 1/4( g_I + g_3 - i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) + 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) - 
			     I * gammaspinmatrix_trace( g0g5 , D ) - 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L1 : // does 1/4( g_I - g_3 + i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) +
			     I * gammaspinmatrix_trace( g0g5 , D ) - 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L2 : // does 1/4( g_I - g_3 - i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) -
			     I * gammaspinmatrix_trace( g0g5 , D ) + 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L3 : // does 1/4( g_I + g_3 + i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) + 
			     gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) +
			     I * gammaspinmatrix_trace( g0g5 , D ) + 
			     I * gammaspinmatrix_trace( g3g0g5 , D ) ) ;
      break ;
    case L4 : // does 1/2 ( g_I + g_3 )
      result[ t ] = 0.5 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) +
			    gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) ) ;
      break ;
    case L5 : // does 1/2 ( g_I - g_3 )
      result[ t ] = 0.5 * ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
			    gammaspinmatrix_trace( GAMMA[ GAMMA_3 ] , D ) ) ;
      break ;
    }
  }

  free( D ) ;
  return result ;
}
