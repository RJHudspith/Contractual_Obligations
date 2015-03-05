/**
   @file momspace_PImunu.c
   @brief computes the momentum-space \f$ \Pi_{\mu\nu}(q) \f$ from our data

   TODO :: projections, cuts, etc etc
 */

#include "common.h"

#include "geometry.h"  // mapping between 0->2Pi BZ that FFTW uses
#include "plan_ffts.h" // for FFTW planning

// enumerated type for the projection
typedef enum {
  TRANS_PROJ ,
  LONG_PROJ ,
  ZERO_PLUS_ONE } projtype ;

// projections
static void
project( struct PIdata *data ,
	 const double **p ,
	 const double *psq ,
	 const projtype PROJ )
{
  const double NORM = 1.0 / (double)( ND - 1 ) ;
  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    const double spsq = ( psq[i] == 0.0 ) ? 1.0 : 1.0 / psq[i] ;
    double sum = 0.0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	const double pmunu = p[i][mu] * p[i][nu] * spsq ;
	double fac = 0.0 ;
	switch( PROJ ) {
	case TRANS_PROJ :
	  fac = ( mu != nu ) ? -pmunu : 1.0 - pmunu ;
	  break ;
	case LONG_PROJ :
	  fac = -pmunu ;
	  break ;
	case ZERO_PLUS_ONE :
	  fac = ( mu == nu ) ? 1.0 - ( ND ) * pmunu : -( ND ) * pmunu ;
	  break ;
	}
	sum += creal( data[i].PI[mu][nu] ) * fac ;
      }
    }
    sum *= ( spsq * NORM ) ;
    printf( "%f %e \n" , psq[i] , sum ) ;
  }
  return ;
}

// \hat{q}_\mu \Pi_{\mu\nu} 
static void
compute_WI( struct PIdata *data ,
	    double **p ,
	    double *psq )
{
  // momenta and correction
  double sum = 0.0 , sum2 = 0.0 ;
  int i ;
  #pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // compute psq and p[i] and compute the correction
    double cache ;
    register double complex epimu = 0.0 ;

    psq[i] = 0.0 ;
    int mu , nu , k[ ND ] = {} ;  // Fourier mode
    TwoPI_mpipi_momconv( k , i , ND ) ;

    for( mu = 0 ; mu < ND ; mu++ ) {
  
      // compute the momenta
      cache = k[ mu ] * Latt.twiddles[ mu ] ;

      // lattice momentum
      p[i][mu] = 2.0 * sin( cache * 0.5 ) ;
      psq[i] += p[i][mu] * p[i][mu] ;

      // WI correction is only for the Vector no?
      epimu = cos( cache * 0.5 ) - I * sin( cache * 0.5 ) ;
      for( nu = 0 ; nu < ND ; nu++ ) {
	// perform multiplicative correction
	data[ i ].PI[mu][nu] *= epimu ;
      }
    }

    register double complex loc_sum = 0.0 , loc_sum2 = 0.0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	// compute ward identities
        loc_sum  += p[i][mu] * data[ i ].PI[mu][nu] ;
	loc_sum2 += data[ i ].PI[mu][nu] * p[i][nu] ;
      }
    }

    sum  = sum  + creal( loc_sum ) * creal( loc_sum ) + \
                  cimag( loc_sum ) * cimag( loc_sum ) ;
    sum2 = sum2 + creal( loc_sum2 ) * creal( loc_sum2 ) + \
                  cimag( loc_sum2 ) * cimag( loc_sum2 ) ;
    // end of i loop
  }

  const double NORM = 1.0 / (double)( LVOLUME * ND ) ;
  printf( "|| p_{mu} Pi_{mu,nu} || :: %e\n" , sum * NORM ) ;
  printf( "|| PI_{mu,nu} p_{nu} || :: %e\n\n" , sum2 * NORM ) ;
}

static void
tmoments( struct PIdata *data ,
	  const projtype PROJ ) 
{
  struct PIdata *ct = malloc( L0 * sizeof( struct PIdata ) ) ;
  struct PIdata *cp = malloc( L0 * sizeof( struct PIdata ) ) ;

  // do the 0 spatial momentum sum first
  int t ;
#pragma omp parallel for private(t)
  for( t = 0 ; t < L0 ; t++ ) {
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	register double complex sum = 0.0 ;
	int i ;
	for( i = 0 ; i < LCU ; i++ ) {
	  sum += data[ i + LCU * t ].PI[ mu ][ nu ] ; ;
	}
	ct[ t ].PI[ mu ][ nu ] = sum ;
      }
    }
  }

  // DFT in time direction, cp is the DFT in the time direction
  const int lt_2 = L0 >> 1 ;
  int pt = 0 ;
#pragma omp parallel for private( pt ) 
  for( pt = -lt_2 ; pt < lt_2 ; pt++ ) {
    double complex twiddles[ L0 ] ;
    const int p_idx = pt + lt_2 ;
    int t ;
    double cache = 0.0 ;
    for( t = 0 ; t < L0 ; t++ ) {
      cache = ( ( pt ) * t * Latt.twiddles[ ND-1 ] ) ;
      twiddles[ t ] = cos( cache ) + I * sin( cache ) ;
    }
    // 
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	register double complex sum = 0.0 ;
	for( t = 0 ; t < L0 ; t++ ) {
	  sum += twiddles[t] * ct[ t ].PI[ mu ][ nu ] ;
	}
	cp[ p_idx ].PI[ mu ][ nu ] = sum ;
      }
    }
    //
  }
  free( ct ) ;

  for( pt = -lt_2 ; pt < lt_2 ; pt++ ) {
    const int pidx = pt + lt_2 ;
    const double psq = 4.0 * sin( 0.5 * pt * Latt.twiddles[ ND - 1 ] ) * sin( 0.5 * pt * Latt.twiddles[ ND - 1 ] ) ;
    const double invpsq = ( psq == 0.0 ) ? 1.0 : 1.0 / psq ;

    double complex sum = 0.0 ;
    int mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      // projections
      switch( PROJ ) {
      case TRANS_PROJ :
	sum += cp[ pidx ].PI[ mu ][ mu ] ;
	break ;
      case LONG_PROJ :
	sum += -cp[ pidx].PI[ ND-1 ][ ND-1 ] ;
	break ;
      case ZERO_PLUS_ONE :
	sum += cp[ pidx ].PI[ mu ][ mu ] - cp[ pt + lt_2 ].PI[ ND-1 ][ ND-1 ] ;
	break ;
      }
    }
    sum *= invpsq / (double)( ND - 1 ) ;

    printf( "%e %e %e \n" , psq , creal( sum ) , cimag( sum ) ) ;
  }

  free( cp ) ;

  return ;
}

void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV )
{
  const projtype PROJ = TRANS_PROJ ;

  // have a look at the ( 0 , 0 , 0 , T ) data TMOMENTS ?
  tmoments( DATA_VV , PROJ ) ;
  //tmoments( DATA_AA , PROJ ) ;

  // if we have FFTW we can unleash it
#ifdef HAVE_FFTW3_H
  double complex *in = malloc( LVOLUME * sizeof( double complex ) ) ;
  double complex *out = malloc( LVOLUME * sizeof( double complex ) ) ;
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  // might be (but I doubt it) slow, most memory efficient way I could think of
  int mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      /////////////////////////////////////////////
      int i ;
      // AA
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	in[ i ] = DATA_AA[ i ].PI[ mu ][ nu ] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	DATA_AA[ i ].PI[ mu ][ nu ] = out[ i ] ;
      }
      // VV
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	in[ i ] = DATA_VV[ i ].PI[ mu ][ nu ] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	DATA_VV[ i ].PI[ mu ][ nu ] = out[ i ] ;
      }
      ///////////////////////////////////////////////
    }
  }
  // free fft memory
  free( in ) ;
  free( out ) ;

  fftw_destroy_plan( forward ) ;
  fftw_destroy_plan( backward ) ;

  fftw_cleanup( ) ;

  // allocate momenta
  double *psq = malloc( LVOLUME * sizeof( double ) ) ;
  double **p = malloc( LVOLUME * sizeof( double* ) ) ;

  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    p[i] = (double*)malloc( ND * sizeof( double ) ) ;
  }

  // check the WI
  compute_WI( DATA_VV , p , psq ) ;

  // projection
  project( DATA_VV , (const double**)p , (const double*)psq , PROJ ) ;

  // free our momenta
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( p[i] ) ;
  }
  free( p ) ;
  free( psq ) ;
#endif

  return ;
}
