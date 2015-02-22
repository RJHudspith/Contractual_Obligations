/**
   @file conserved_local.c
   @brief conserved-local Wilson currents
 */

#include "common.h"

#include "plan_ffts.h"

// struct to store PI_{\mu\nu}
struct fftdata {
  double complex PI[ ND ][ ND ] ;
} ;

static double complex
contract_conserved_local( struct fftdata *FFTDATA_VV ,
			  struct fftdata *FFTDATA_AA ,
			  struct site *lat ,
			  const struct spinor S1 ,
			  const struct spinor S1UP ) 
{
  return 0.0 ;
}

int
conserved_local( FILE *fprop1 ,
		 const int header )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS ) ;

  // over the whole volume is not as expensive as you may think
  struct fftdata *FFTDATA_VV = malloc( LVOLUME * sizeof( struct fftdata ) ) ;
  struct fftdata *FFTDATA_AA = malloc( LVOLUME * sizeof( struct fftdata ) ) ;
  
#ifdef HAVE_FFTW3_H
  double complex *in = malloc( LVOLUME * sizeof( double complex ) ) ;
  double complex *out = malloc( LVOLUME * sizeof( double complex ) ) ;
  fftw_plan forwardVV , backward ;
  small_create_plans_DFT( forward , backward , in , out , ND ) ;
#endif

  // and our spinor
  struct spinor *S1 = malloc( VOL3 * sizeof( struct spinor ) ) ;
  struct spinor *S1UP = malloc( VOL3 * sizeof( struct spinor ) ) ;

#ifdef HAVE_FFTW3_H
  // might be (but I doubt it) slow, most memory efficient way I could think of
  int mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      // VV
      int i ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	in[ i ] = FFTDATA_VV[ i ].PI[ mu ][ nu ] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	FFTDATA_VV[ i ].PI[ mu ][ nu ] = out[ i ] ;
      }
      // AA
      int i ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	in[ i ] = FFTDATA_AA[ i ].PI[ mu ][ nu ] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	FFTDATA_AA[ i ].PI[ mu ][ nu ] = out[ i ] ;
      }
    }
  }

  free( in ) ;
  free( out ) ;

  fftw_destroy_plan( forward ) ;
  fftw_destroy_plan( backward ) ;

  fftw_cleanup( ) ;
#endif

  // now, the CL data will be in the 0 -> 2\Pi BZ so we have to be careful
  // in geometry there are some routines we can use to translate between this
  // and the -Pi -> Pi BZ which we usually live in

  // do a projection to get transverse, longitudinal & their sum and also sum cuts

  // free the CL data
  free( FFTDATA_VV ) ;
  free( FFTDATA_AA ) ;

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S2 ) ;

  return SUCCESS ;
}
