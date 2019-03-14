/**
   @file bar_projections.c
   @brief spin-1/2 and spin-3/2 projections
 */
#include "common.h"

#include "bar_projections.h" // alphabetising
#include "gammas.h"          // gamma_mmul()
#include "spinmatrix_ops.h"  // matrix operations

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
  if( corr_malloc( (void**)&tmp1 , ALIGNMENT ,
		   NSNS*sizeof( double complex ) ) != 0 ||
      corr_malloc( (void**)&tmp2 , ALIGNMENT ,
		   NSNS*sizeof( double complex ) ) != 0 ) {
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
  if( corr_malloc( (void**)&tmp1 , ALIGNMENT ,
		   NSNS*sizeof( double complex ) ) != 0 ||
      corr_malloc( (void**)&tmp2 , ALIGNMENT ,
		   NSNS*sizeof( double complex ) ) != 0 ) {
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

// select matrix G_{ik} with dirac indices i (row) k (cols)
static void
set_Gik( double complex *Gik , 
	 const struct mcorr **corr ,
	 const size_t i ,
	 const size_t k ,
	 const size_t pidx ,
	 const size_t t ,
	 const size_t NGAMS )
{
  size_t d ;
#if NS == 4
  for( d = 0 ; d < NS ; d++ ) {
    Gik[ 0 + d * NS ] = corr[ k + i*NGAMS ][ 0 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 1 + d * NS ] = corr[ k + i*NGAMS ][ 1 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 2 + d * NS ] = corr[ k + i*NGAMS ][ 2 + d * NS ].mom[ pidx ].C[ t ] ;
    Gik[ 3 + d * NS ] = corr[ k + i*NGAMS ][ 3 + d * NS ].mom[ pidx ].C[ t ] ;
  }
#else
  for( d = 0 ; d < NSNS ; d++ ) {
    Gik[ d ] = corr[ k + i*NGAMS ][ d ].mom[ pidx ].C[ t ] ;
  }
#endif
}

// check if a gamma index can be spin-projected
static int
check_idx( const size_t idx )
{
  switch( (int)idx%16 ) {
  case GAMMA_X : case AX : case TXT :
    return 0 ;
  case GAMMA_Y : case AY : case TYT :
    return 1 ;
  case GAMMA_Z : case AZ : case TZT :
    return 2 ;
  case GAMMA_T : case IDENTITY : case GAMMA_5 :
  case AT :
    // not actually sure about the tensors
  case TXY : case TYZ : case TZX : 
    return GAMMA_T ;
  }
  return GAMMA_T ;
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
	      const spinhalf projection ,
	      const size_t NGAMS ) 
{
  // spin projections are for spatial indices only!!
  if( check_idx(i) == GAMMA_T || check_idx(k) == GAMMA_T ) {
    set_Gik( Gik , corr , i , k , p , t , NGAMS ) ;
    return ;
  }
  
  // do the spin projections
  switch( projection ) {
  case OneHalf_11 :
    spinproj( Gik , P11 , corr , i , k , t , GAMMA , momentum , p , NGAMS ) ;
    break ;
  case OneHalf_12 :
    spinproj( Gik , P12 , corr , i , k , t , GAMMA , momentum , p , NGAMS ) ;
    break ;
  case OneHalf_21 :
    spinproj( Gik , P21 , corr , i , k , t , GAMMA , momentum , p , NGAMS ) ;
    break ;
  case OneHalf_22 :
    spinproj( Gik , P22 , corr , i , k , t , GAMMA , momentum , p , NGAMS ) ;
    break ;
  case ThreeHalf :
    spinproj( Gik , P32 , corr , i , k , t , GAMMA , momentum , p , NGAMS ) ;
    break ;
  case NONE :
    set_Gik( Gik , corr , i , k , p , t , NGAMS ) ;
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
		const spinhalf spin_proj ,
		const size_t NGAMS ) 
{
  // D is a temporary spinmatrix
  double complex *D = NULL , *result = NULL ;

  // allocate temporaries
  corr_malloc( (void**)&D , ALIGNMENT , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&result , ALIGNMENT , LT * sizeof( double complex ) ) ;

  // build these just in case
  struct gamma gxg5 , gtgxg5 , gyg5 , gtgyg5 , gzg5 , gtgzg5 ;
  // x-projection
  gamma_mmul( &gxg5 , GAMMA[ GAMMA_Z ] , GAMMA[ GAMMA_5 ] ) ;
  gamma_mmul( &gtgxg5 , GAMMA[ GAMMA_T ] , gxg5 ) ;
  // y-projection
  gamma_mmul( &gyg5 , GAMMA[ GAMMA_Z ] , GAMMA[ GAMMA_5 ] ) ;
  gamma_mmul( &gtgyg5 , GAMMA[ GAMMA_T ] , gyg5 ) ;
  // z-projection
  gamma_mmul( &gzg5 , GAMMA[ GAMMA_Z ] , GAMMA[ GAMMA_5 ] ) ;
  gamma_mmul( &gtgzg5 , GAMMA[ GAMMA_T ] , gzg5 ) ;
  
  // loop times 
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {

    // set open dirac matrix doing the spin projection if desired
    spin_project( D , corr , GAMMA , momentum , 
		  GSRC , GSNK , p , t , spin_proj , NGAMS ) ; 
    
    // Parity projections are ::
    // L0_i = ( 1 + \gamma_t ) * ( 1 + i \gamma_i \gamma_5 )/4
    // L1_i = ( 1 + \gamma_t ) * ( 1 - i \gamma_i \gamma_5 )/4
    // L2_i = ( 1 - \gamma_t ) * ( 1 + i \gamma_i \gamma_5 )/4
    // L3_i = ( 1 - \gamma_t ) * ( 1 - i \gamma_i \gamma_5 )/4
    // Note that L^{i}_j.L^{i}_j = L^{i}_j
    // BUT
    // L_4.L5 = 0, L_{0/1}.L_{1/0,2,3,5} = 0, L_{2/3}.L_{0,1,3/2,4} = 0
    switch( parity_proj ) {
    case L0 :
      result[ t ] = ( +3*gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D )
		      +3*gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D )  
		      +I * gammaspinmatrix_trace( gxg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgxg5 , D ) 
		      +I * gammaspinmatrix_trace( gyg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgyg5 , D ) 
		      +I * gammaspinmatrix_trace( gzg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgzg5 , D ) ) / 12. ;
      break ;
    case L1 :
      result[ t ] = ( +3*gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D )
		      +3*gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D )  
		      -I * gammaspinmatrix_trace( gxg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgxg5 , D ) 
		      -I * gammaspinmatrix_trace( gyg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgyg5 , D ) 
		      -I * gammaspinmatrix_trace( gzg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgzg5 , D ) ) / 12. ;
      break ;
    case L2 :
      result[ t ] = ( +3*gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D )
		      -3*gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D )
		      +I * gammaspinmatrix_trace( gxg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgxg5 , D ) 
		      +I * gammaspinmatrix_trace( gyg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgyg5 , D )
		      +I * gammaspinmatrix_trace( gzg5 , D ) 
		      -I * gammaspinmatrix_trace( gtgzg5 , D ) ) / 12. ;
      break ;
    case L3 :
      result[ t ] = ( +3*gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D )
		      -3*gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D )
		      -I * gammaspinmatrix_trace( gxg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgxg5 , D ) 
		      -I * gammaspinmatrix_trace( gyg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgyg5 , D ) 
		      -I * gammaspinmatrix_trace( gzg5 , D ) 
		      +I * gammaspinmatrix_trace( gtgzg5 , D ) ) / 12. ;
      break ;
    case L4 : // does 1/2 ( g_I + g_t )
      result[ t ] = ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) +
		      gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) ) / 2. ;
      break ;
    case L5 : // does 1/2 ( g_I - g_t )
      result[ t ] = ( gammaspinmatrix_trace( GAMMA[ IDENTITY ] , D ) - 
		      gammaspinmatrix_trace( GAMMA[ GAMMA_T ] , D ) ) / 2. ;
      break ;
    }
  }

  free( D ) ;
  return result ;
}

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
    if( corr_malloc( (void**)&left , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 ||
	corr_malloc( (void**)&right , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 || 
	corr_malloc( (void**)&pslash , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0  ) {
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
    if( corr_malloc( (void**)&pslash , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 ) {
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
    if( corr_malloc( (void**)&pslash , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 ) {
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
    if( corr_malloc( (void**)&left , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 ||
	corr_malloc( (void**)&right , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0 || 
	corr_malloc( (void**)&pslash , ALIGNMENT , NSNS*sizeof( double complex ) ) != 0  ) {
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
	  const size_t pidx ,
	  const size_t NGAMS )
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
  if( corr_malloc( (void**)&pij , ALIGNMENT ,
		   NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }
  if( corr_malloc( (void**)&tmp , ALIGNMENT ,
		   NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }
  if( corr_malloc( (void**)&Gjk , ALIGNMENT ,
		   NSNS * sizeof( double complex ) ) != 0 ) {
    flag = FAILURE ;
    goto memfree ;
  }

  // perform contraction
  const size_t gi = check_idx( i ) ;
  for( j = 0 ; j < ND ; j++ ) {
    set_Gik( Gjk , corr , j , k , pidx , t , NGAMS ) ;
    p( pij , gi , j , GAMMA , momentum , pidx ) ;
    spinmatrix_multiply( tmp , pij , Gjk ) ;
    atomic_add_spinmatrices( Gik , tmp ) ;
  }

 memfree :
    
  free( tmp ) ;
  free( pij ) ;
  free( Gjk ) ;

  return flag ;
}
