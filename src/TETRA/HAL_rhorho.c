/**
   @file dibaryon_contractions.c
   @brief SU(2) dibaryon contraction code
 */
#include "common.h"

#include "gammas.h"
#include "geometry.h"
#include "Ospinor.h"
#include "spinmatrix_ops.h" 
#include "spinor_ops.h"     // transpose_spinor()

// does Tr[ S1T G1 S2 G2 ( S3T G3 S4 G4 )^\dagger ]
static inline double complex
single_tr( struct spinmatrix block1 ,
	   struct spinmatrix block2 )
{
  return trace_prod_spinmatrices( (const void*)block1.D ,
				  (const void*)block2.D ) ;
}

// does Tr[ S1T G1 S2 G2 ] x Tr[ ( S3T G3 S4 G4 )^\dagger ]
static inline double complex
double_tr( struct spinmatrix block1 ,
	   struct spinmatrix block2 )
{
  return
    spinmatrix_trace( (const void*)block1.D ) *
    ( spinmatrix_trace( (const void*)block2.D ) ) ;
}

// from lexicographical index abcd get index a,b,c,d --> only correct for NC=2
static inline void
get_abcd( size_t *a , size_t *b , size_t *ap , size_t *bp , const size_t abcd )
{
  *a  = abcd >> 3 ;
  *b  = ( abcd >> 2 )&1 ;
  *ap = ( abcd >> 1 )&1 ;
  *bp = abcd&1 ;
  return ;
}

// given a,b,ap,bp return lexicographical index abcd
static inline size_t
get_idx( const size_t a , const size_t b , const size_t ap , const size_t bp )
{
  return bp + NC*( ap + NC*( b + NC*a ) ) ;
}

// do all the precomputations
static void
precompute_blocks( struct spinmatrix **blk11 ,
		   struct spinmatrix **blk12 ,
		   struct spinmatrix **blk21 ,
		   struct spinmatrix **blk22 ,
		   const struct spinor *S ,
		   const struct gamma *GAMMAS ,
		   const size_t GSRC ,
		   const size_t GSNK )
{
  // gamma precomputations
  const struct gamma G1   = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
  const struct gamma tG1t = gt_Gdag_gt( G1 , GAMMAS[ GAMMA_T ] ) ;

  // gamma precomputations
  const struct gamma G2   = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
  const struct gamma tG2t = gt_Gdag_gt( G2 , GAMMAS[ GAMMA_T ] ) ;
  
  size_t i ;
  // contract over the volume
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
  
    // precomputations -> swap color and dirac indices to expose spinmatrices
    struct Ospinor OST = spinor_to_Ospinor( transpose_spinor( S[i] ) ) ;
    struct Ospinor OS  = spinor_to_Ospinor( S[i] ) ;
  
    // multiply left and right by gammas into tmp with color indices ab
    // and multiply on the left by the transposed prop and save result
    struct spinmatrix tmp ;
    size_t a , b , ap , bp , abcd ;
    for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
      get_abcd( &a , &b , &ap , &bp , abcd ) ;
      // g1g1 case
      tmp = OS.C[a][b] ;
      gamma_spinmatrix_lr( &tmp , G1 , tG1t ) ;
      spinmatrix_multiply( (void*)blk11[i][abcd].D ,
			   (void*)OST.C[ap][bp].D , (void*)tmp.D ) ;
      // g1g2 case
      tmp = OS.C[a][b] ;
      gamma_spinmatrix_lr( &tmp , G1 , tG2t ) ;
      spinmatrix_multiply( (void*)blk12[i][abcd].D ,
			   (void*)OST.C[ap][bp].D , (void*)tmp.D ) ;
      // g2g1 case
      tmp = OS.C[a][b] ;
      gamma_spinmatrix_lr( &tmp , G2 , tG1t ) ;
      spinmatrix_multiply( (void*)blk21[i][abcd].D ,
			   (void*)OST.C[ap][bp].D , (void*)tmp.D ) ;
      // g2g2 case
      tmp = OS.C[a][b] ;
      gamma_spinmatrix_lr( &tmp , G2 , tG2t ) ;
      spinmatrix_multiply( (void*)blk22[i][abcd].D ,
			   (void*)OST.C[ap][bp].D , (void*)tmp.D ) ;
    }
  }
  return ;
}

// perform the required FFTs
static void
FFT_blocks( double complex *in ,
	    double complex *out ,
	    struct spinmatrix **blk11 ,
	    struct spinmatrix **blk12 ,
	    struct spinmatrix **blk21 ,
	    struct spinmatrix **blk22 ,
	    const fftw_plan forward ,
	    const fftw_plan backward )
{
  // TODO :: this really should be threaded or something
#pragma omp single
  {
    size_t abcd , i , d ;
    for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
      // loop dirac indices
      for( d = 0 ; d < NSNS ; d++ ) {
	const size_t d1 = d/NS ;
	const size_t d2 = d%NS ;
	// set "in"
	for( i = 0 ; i < LCU ; i++ ) {
	  in[ i ] = blk11[ i ][ abcd ].D[ d1 ][ d2 ] ;
	}
	fftw_execute( forward ) ;
	for( i = 0 ; i < LCU ; i++ ) {
	  blk11[ i ][ abcd ].D[ d1 ][ d2 ] = out[i] ;
	  in[ i ] = blk12[ i ][ abcd ].D[ d1 ][ d2 ] ;
	}
	fftw_execute( forward ) ;
	for( i = 0 ; i < LCU ; i++ ) {
	  blk12[ i ][ abcd ].D[ d1 ][ d2 ] = out[i] ;
	  out[ i ] = blk21[ i ][ abcd ].D[ d1 ][ d2 ] ;
	}
	// these two need to be backward transforms as the
	// convolution wants a "-p" transform
	fftw_execute( backward ) ;
	for( i = 0 ; i < LCU ; i++ ) {
	  blk21[ i ][ abcd ].D[ d1 ][ d2 ] = in[i] ;
	  out[ i ] = blk22[ i ][ abcd ].D[ d1 ][ d2 ] ;
	}
	fftw_execute( backward ) ;
	for( i = 0 ; i < LCU ; i++ ) {
	  blk22[ i ][ abcd ].D[ d1 ][ d2 ] = in[i] ;
	}
      }
    }
  }
  return ;
}

// do the contractions
static void
contract_diagrams( double complex *out ,
		   const struct spinmatrix **blk11 ,
		   const struct spinmatrix **blk12 ,
		   const struct spinmatrix **blk21 ,
		   const struct spinmatrix **blk22 ,
		   const int r[ ND ] )
{
  // do all of the 4! contractions summing over the 16 color combinations
  // note that I use the identity that G^T for gamma_i is G
  // and then convolve
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // compute the shifted spacing
    const size_t sep = compute_spacing( r , i , ND-1 ) ;
    
    // point to the blocks
    const struct spinmatrix *b11 = blk11[sep] ;
    const struct spinmatrix *b12 = blk12[sep] ;
    const struct spinmatrix *b21 = blk21[i] ;
    const struct spinmatrix *b22 = blk22[i] ;
    // summation of all the diagrams goes into sum
    register double complex sum = 0.0 ;
    size_t abcd , a , b , ap , bp ;
    for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
      get_abcd( &a , &b , &ap , &bp , abcd ) ;
      // diagram 1/24
      sum += double_tr( b11[get_idx(a,ap,b,bp)] , b22[get_idx(a,ap,b,bp)] ) ;
      // diagram 2/24
      sum -= double_tr( b11[get_idx(b,ap,a,bp)] , b22[get_idx(a,ap,b,bp)] ) ;
      // diagram 3/24
      sum += single_tr( b11[get_idx(b,ap,a,ap)] , b22[get_idx(a,bp,b,bp)] ) ;
      // diagram 4/24
      sum -= single_tr( b11[get_idx(a,ap,b,ap)] , b22[get_idx(a,bp,b,bp)] ) ;
      // diagram 5/24
      sum -= single_tr( b11[get_idx(b,ap,a,bp)] , b22[get_idx(a,ap,b,bp)] ) ;
      // diagram 6/24
      sum += single_tr( b11[get_idx(a,ap,b,bp)] , b22[get_idx(a,ap,b,bp)] ) ;
      // diagram 7/24
      sum -= double_tr( b11[get_idx(a,ap,b,bp)] , b22[get_idx(b,ap,a,bp)] ) ;
      // diagram 8/24
      sum += double_tr( b11[get_idx(b,ap,a,bp)] , b22[get_idx(b,ap,a,bp)] ) ;
      // diagram 9/24
      sum -= single_tr( b12[get_idx(a,ap,b,bp)] , b21[get_idx(b,ap,a,bp)] ) ;
      // diagram 10/24
      sum += single_tr( b12[get_idx(b,ap,a,bp)] , b21[get_idx(b,ap,a,bp)] ) ;
      // diagram 11/24
      sum += single_tr( b11[get_idx(b,bp,a,bp)] , b22[get_idx(a,ap,b,ap)] ) ;
      // diagram 12/24
      sum -= single_tr( b11[get_idx(a,bp,b,bp)] , b22[get_idx(a,ap,b,ap)] ) ;
      //diagram 13/24
      sum -= single_tr( b12[get_idx(a,ap,b,ap)] , b21[get_idx(a,bp,b,bp)] ) ;
      // diagram 14/24
      sum += single_tr( b11[get_idx(a,ap,b,ap)] , b22[get_idx(b,bp,a,bp)] ) ;
      // diagram 15/24
      sum += single_tr( b12[get_idx(a,ap,b,bp)] , b21[get_idx(a,ap,b,bp)] ) ;
      // diagram 16/24
      sum -= single_tr( b12[get_idx(b,ap,a,bp)] , b21[get_idx(a,ap,b,bp)] ) ;
      // diagram 17/24
      sum += double_tr( b12[get_idx(a,ap,b,bp)] , b21[get_idx(a,ap,b,bp)] ) ;
      // diagram 18/24
      sum -= double_tr( b12[get_idx(b,ap,a,bp)] , b21[get_idx(a,ap,b,bp)] ) ;
      // diagram 19/24
      sum += single_tr( b11[get_idx(b,ap,a,bp)] , b22[get_idx(b,ap,a,bp)] ) ;
      // diagram 20/24
      sum -= single_tr( b11[get_idx(a,ap,b,bp)] , b22[get_idx(b,ap,a,bp)] ) ;
      // diagram 21/24
      sum -= single_tr( b12[get_idx(a,bp,b,bp)] , b21[get_idx(a,ap,b,ap)] ) ;
      // diagram 22/24
      sum += single_tr( b12[get_idx(b,bp,a,bp)] , b21[get_idx(a,ap,b,ap)] ) ;
      // diagram 23/24
      sum -= double_tr( b12[get_idx(a,ap,b,bp)] , b21[get_idx(b,ap,a,bp)] ) ;
      // diagram 24/24
      sum += double_tr( b12[get_idx(b,ap,a,bp)] , b21[get_idx(b,ap,a,bp)] ) ;
    }
    out[i] = sum ;
  }
  return ;
}

// perform the HAL-QCD-style rho-rho scattering contraction
int
HALrhorho_contract( double complex *in ,
		    double complex *out ,
		    fftw_plan forward ,
		    fftw_plan backward ,
		    struct spinmatrix **blk11 ,
		    struct spinmatrix **blk12 ,
		    struct spinmatrix **blk21 ,
		    struct spinmatrix **blk22 ,
		    const struct spinor *S ,
		    const struct gamma *GAMMAS ,
		    const size_t GSRC ,
		    const size_t GSNK ,
		    const int nmom[1] ,
		    const struct veclist *list )
{  
  // precompute all the required blocks
  precompute_blocks( blk11 , blk12 , blk21 , blk22 , S , GAMMAS , GSRC , GSNK ) ;

#ifdef HAVE_FFTW3_H
  // perform the ffts
  FFT_blocks( in , out , blk11 , blk12 , blk21 , blk22 , forward , backward ) ;
  
  // perform the contractions
  int dummy[ ND ] = { 0 , 0 , 0 , 0 } ;
  contract_diagrams( out ,
		     (const struct spinmatrix**)blk11 ,
		     (const struct spinmatrix**)blk12 ,
		     (const struct spinmatrix**)blk21 ,
		     (const struct spinmatrix**)blk22 ,
		     dummy ) ;
  
  // fft back
  #pragma omp single
  {
    fftw_execute( backward ) ;
  }

  // take the FFT norm
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    in[i] /= (double)LCU ;
  }

  // ok in here we actually want both the second-order derivative
  // grad.\phi(r) approximated by f(x/y/z+1)+f(x/y/z-1)-2f(x/y/z) at each "r" vector point
  // into "out" -> "in" still contains \phi(r)
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    const size_t s1 = gen_shift( i , 0 ) , s2 = gen_shift( i , -1 ) ;
    out[i] = in[ s1 ] + in[ s2 ] - 2*in[ i ] ;
    const size_t s3 = gen_shift( i , 1 ) , s4 = gen_shift( i , -2 ) ;
    out[i] += in[ s3 ] + in[ s4 ] - 2*in[ i ] ;
    const size_t s5 = gen_shift( i , 2 ) , s6 = gen_shift( i , -3 ) ;
    out[i] += in[ s5 ] + in[ s6 ] - 2*in[ i ] ;
  }
  
#else
  size_t r ;
  for( r = 0 ; r < nmom[0] ; r++ ) {
    const int cast[ ND ] = { (int)list[r].MOM[0] , (int)list[r].MOM[1] ,
			     (int)list[r].MOM[2] , (int)list[r].MOM[3] } ;
    contract_diagrams( out ,
		       (const struct spinmatrix**)blk11 ,
		       (const struct spinmatrix**)blk12 ,
		       (const struct spinmatrix**)blk21 ,
		       (const struct spinmatrix**)blk22 ,
		       cast ) ;
    // sum all of "in" into "out array
    double complex sum = 0 ;
    size_t site ;
    #pragma omp single
    {
      for( site = 0 ; site < LCU ; site++ ) {
	sum = sum + out[ site ] ;
      }
      in[ list[ r ].idx ] = sum ;
    }
  }
#endif
  
  return SUCCESS ;
}
