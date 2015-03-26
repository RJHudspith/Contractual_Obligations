/**
   @file momspace_PImunu.c
   @brief computes the momentum-space \f$ \Pi_{\mu\nu}(q) \f$ from our data
 */

#include "common.h"

#include "cut_routines.h"       // cuts
#include "plan_ffts.h"          // for FFTW planning
#include "PImunu_projections.h" // projection codes
#include "WardIdentity.h"       // compute_psq

// FFT our mu,nu data
#ifdef HAVE_FFTW3_H
static void
FFT_PImunu( struct PIdata *DATA ,
	    double complex *in , 
	    double complex *out ,
	    const fftw_plan forward ,
	    const fftw_plan backward ) 
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
#endif

//
void
momspace_PImunu( struct PIdata *AA ,
		 struct PIdata *VV ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ,
		 const current_type current )
{
  // if we have FFTW we can unleash it
#ifdef HAVE_FFTW3_H
  // temporary space
  double complex *in = malloc( LVOLUME * sizeof( double complex ) ) ;
  double complex *out = malloc( LVOLUME * sizeof( double complex ) ) ;

  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  // might be (but I doubt it) slow, most memory efficient way I could think of
  FFT_PImunu( AA , in , out , forward , backward ) ;
  FFT_PImunu( VV , in , out , forward , backward ) ;

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

  momspace_data( AA , (const double **)p , psq , list , 
		 NMOM , outfile , current , AXIAL ) ;

  momspace_data( VV , (const double **)p , psq , list , 
		 NMOM , outfile , current , VECTOR ) ;

  // free our momenta
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    free( p[i] ) ;
  }
  free( p ) ;
  free( psq ) ;

  // free the momentum list
  free( NMOM ) ;
  free( (void*)list ) ;

#else
  printf( "NON-fftw routines not supported yet \n" ) ;
#endif
  return ;
}
