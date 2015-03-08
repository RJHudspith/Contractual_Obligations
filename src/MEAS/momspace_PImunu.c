/**
   @file momspace_PImunu.c
   @brief computes the momentum-space \f$ \Pi_{\mu\nu}(q) \f$ from our data

   TODO :: cuts, etc etc
 */

#include "common.h"

#include "cut_routines.h" // cuts
#include "geometry.h"     // mapping between 0->2Pi BZ that FFTW uses
#include "plan_ffts.h"    // for FFTW planning
#include "WardIdentity.h" // compute the ward identity

// projections
static void
project( double *trans ,
	 double *longitudinal ,
	 const struct PIdata *data ,
	 const double **p ,
	 const double *psq ,
	 const struct veclist *list ,
	 const int NMOM )
{
  const double NORM = 1.0 / (double)( ND - 1 ) ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {
    const int list_idx = list[i].idx ;
    const double spsq = ( psq[i] == 0.0 ) ? 1.0 : 1.0 / psq[i] ;
    double sumtrans = 0.0 , sumlong = 0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	const double pmunu = p[i][mu] * p[i][nu] * spsq ;
	double fac = ( mu != nu ) ? -pmunu : 1.0 - pmunu ;
	sumtrans += creal( data[list_idx].PI[mu][nu] ) * fac ;
	sumlong  += -creal( data[list_idx].PI[mu][nu] ) * pmunu ;
      }
    }
    trans[ i ] = sumtrans * ( spsq * NORM ) ;
    longitudinal[ i ] = sumlong * spsq ;
  }
  return ;
}

// FFT our data
static void
FFT_PImunu( struct PIdata *DATA ,
	    double complex *in , 
	    double complex *out ,
	    fftw_plan forward ,
	    fftw_plan backward ) 
{
  int mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      int i ;
      // AA
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	in[ i ] = DATA[ i ].PI[ mu ][ nu ] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	DATA[ i ].PI[ mu ][ nu ] = out[ i ] ;
      }
    }
  }
  return ;
}

//
void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV ,
		 const struct cut_info CUTINFO )
{
  // if we have FFTW we can unleash it
#ifdef HAVE_FFTW3_H
  double complex *in = malloc( LVOLUME * sizeof( double complex ) ) ;
  double complex *out = malloc( LVOLUME * sizeof( double complex ) ) ;
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  // might be (but I doubt it) slow, most memory efficient way I could think of
  FFT_PImunu( DATA_AA , in , out , forward , backward ) ;
  FFT_PImunu( DATA_VV , in , out , forward , backward ) ;

  // free fft memory
  free( in ) ;
  free( out ) ;

  fftw_destroy_plan( forward ) ;
  fftw_destroy_plan( backward ) ;

  fftw_cleanup( ) ;

  // here is where the call to the cuts routine goes
  int *NMOM = malloc( sizeof( int ) ) ;
  const struct veclist *list = compute_veclist( NMOM , CUTINFO ,
						ND , GLU_FALSE ) ;

  // allocate momenta
  double *psq = malloc( NMOM[0] * sizeof( double ) ) ;
  double **p = malloc( NMOM[0] * sizeof( double* ) ) ;

  int i ;
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    p[i] = (double*)malloc( ND * sizeof( double ) ) ;
  }

  // precompute momenta
  compute_p_psq( p , psq , list , NMOM[0] ) ;

  // correct the WI on the vector
  correct_WI( DATA_VV , CORR_MU , list , NMOM[0] ) ;
  compute_WI( DATA_VV , (const double**)p , list , NMOM[0] ) ;

  double *trans = malloc( NMOM[0] * sizeof( double ) ) ;
  double *longitudinal = malloc( NMOM[0] * sizeof( double ) ) ;
 
  // projection
  project( trans , longitudinal , DATA_VV , 
	   (const double**)p , (const double*)psq ,
	   list , NMOM[0] ) ;

  
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    printf( "%f %f \n" , psq[i] , trans[i] ) ;
  }

  free( NMOM ) ;

  // free our transverse and longitudinal components
  free( trans ) ;
  free( longitudinal ) ;

  // free our momenta
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    free( p[i] ) ;
  }
  free( p ) ;
  free( psq ) ;

#else
  printf( "NON-fftw routines not supported yet \n" ) ;
#endif
  return ;
}
