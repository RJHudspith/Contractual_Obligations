/**
   @file momspace_PImunu.c
   @brief computes the momentum-space \f$ \Pi_{\mu\nu}(q) \f$ from our data
 */

#include "common.h"

#include "cut_output.h"   // cut output format
#include "cut_routines.h" // cuts
#include "geometry.h"     // mapping between 0->2Pi BZ that FFTW uses
#include "plan_ffts.h"    // for FFTW planning
#include "WardIdentity.h" // compute the ward identity

// projections should go in their own file right?
static void
project( const struct PIdata *data ,
	 const double **p ,
	 const double *psq ,
	 const struct veclist *list ,
	 const int *__restrict NMOM ,
	 const char *outfile )
{
  double *trans = malloc( NMOM[0] * sizeof( double ) ) ;
  double *longitudinal = malloc( NMOM[0] * sizeof( double ) ) ;

  const double NORM = 1.0 / (double)( ND - 1 ) ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM[0] ; i++ ) {
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

  // write out a file
  char str[ 256 ] ;
  sprintf( str , "%s.trans" , outfile ) ;
  write_momspace_data( str , NMOM , trans , list , ND ) ;

  sprintf( str , "%s.long" , outfile ) ;
  write_momspace_data( str , NMOM , longitudinal , list , ND ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    trans[ i ] += longitudinal[ i ] ; // put the 0+1 contribution in too
  }

  sprintf( str , "%s.transPlong" , outfile ) ;
  write_momspace_data( str , NMOM , longitudinal , list , ND ) ;

  // free the projected data
  free( trans ) ;
  free( longitudinal ) ;

  return ;
}

// we have the zero at index [0] as the FFT is in the 0->2Pi BZ
static void
subtract_zeromom( struct PIdata *__restrict data )
{
  // precompute subtraction point
  double complex sub[ ND ][ ND ] = {} ;
  int i , j ;
  for( i = 0 ; i < ND ; i++ ) {
    for( j = 0 ; j < ND ; j++ ) {
      sub[ i ][ j ] = data[ 0 ].PI[ i ][ j ] ;
    }
  }
  // and subtract the noisy zero
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	data[ i ].PI[ mu ][ nu ] -= sub[ mu ][ nu ] ;
      }
    }
  }
  return ;
}

// FFT our mu,nu data
#ifdef HAVE_FFTW3_H
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
#endif

//
void
momspace_PImunu( struct PIdata *DATA_AA ,
		 struct PIdata *DATA_VV ,
		 const struct cut_info CUTINFO ,
		 const char *outfile )
{
  // if we have FFTW we can unleash it
#ifdef HAVE_FFTW3_H
  // temporary space
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
 
  // projection
  project( DATA_VV , (const double**)p , (const double*)psq ,
	   list , NMOM , outfile ) ;

  // free the momentum list
  free( NMOM ) ;
  free( (void*)list ) ;

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
