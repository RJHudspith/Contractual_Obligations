/**
   @file momspace_PImunu.c
   @brief computes the momentum-space \f$ \Pi_{\mu\nu}(q) \f$ from our data

   TODO :: cuts, etc etc
 */

#include "common.h"

#include "geometry.h"     // mapping between 0->2Pi BZ that FFTW uses
#include "plan_ffts.h"    // for FFTW planning
#include "WardIdentity.h" // compute the ward identity

// projections
static void
project( double *trans ,
	 double *longitudinal ,
	 const struct PIdata *data ,
	 const double **p ,
	 const double *psq )
{
  const double NORM = 1.0 / (double)( ND - 1 ) ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    const double spsq = ( psq[i] == 0.0 ) ? 1.0 : 1.0 / psq[i] ;
    double sumtrans = 0.0 , sumlong = 0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	const double pmunu = p[i][mu] * p[i][nu] * spsq ;
	double fac = ( mu != nu ) ? -pmunu : 1.0 - pmunu ;
	sumtrans += creal( data[i].PI[mu][nu] ) * fac ;
	sumlong  += -creal( data[i].PI[mu][nu] ) * pmunu ;
      }
    }
    trans[ i ] = sumtrans * ( spsq * NORM ) ;
    longitudinal[ i ] = sumlong * spsq ;
  }
  return ;
}

//
void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV )
{
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

  // here is where the call to the cuts routine goes

  // allocate momenta
  double *psq = malloc( LVOLUME * sizeof( double ) ) ;
  double **p = malloc( LVOLUME * sizeof( double* ) ) ;

  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    p[i] = (double*)malloc( ND * sizeof( double ) ) ;
  }

  // precompute momenta
  compute_p_psq( p , psq ) ;

  // correct the WI on the vector
  correct_WI( DATA_VV , CORR_MU ) ;
  compute_WI( DATA_VV , (const double**)p ) ;

  double *trans = malloc( LVOLUME * sizeof( double ) ) ;
  double *longitudinal = malloc( LVOLUME * sizeof( double ) ) ;
 
  // projection
  project( trans , longitudinal , DATA_VV , 
	   (const double**)p , (const double*)psq ) ;

  free( trans ) ;
  free( longitudinal ) ;

  // free our momenta
  for( i = 0 ; i < LVOLUME ; i++ ) {
    free( p[i] ) ;
  }
  free( p ) ;
  free( psq ) ;

#else
  printf( "NON-fftw routines not supported yet \n" ) ;
#endif
  return ;
}
