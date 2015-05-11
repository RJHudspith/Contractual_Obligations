/**
   @file tmoments_PImunu.c
   @brief computation of \f$PI(q^2)\f$ where \f$q=(0,0,0,q_t)\f$

   TODO :: readers & writers for this and momspace PImunu data
 */

#include "common.h"

#include "correlators.h"        // write out the correlator file
#include "cut_routines.h"       // zero_veclist
#include "WardIdentity.h"       // psq calculation
#include "PImunu_projections.h" // for the momentum space projections

// DFT in time direction
// returns veclist so we can reuse the IO code DFT in the 0->2Pi BZ
static struct veclist*
DFT( struct PIdata *cpAA ,
     struct PIdata *cpVV ,
     const struct mcorr **ctAA ,
     const struct mcorr **ctVV )
{
  struct veclist *list = malloc( L0 * sizeof( struct veclist ) ) ;

  const int lt_2 = L0 >> 1 ;
  int pt = 0 ;
#pragma omp parallel for private( pt ) 
  for( pt = 0 ; pt < L0 ; pt++ ) {
 
    // set up the veclist
    int mu , nu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      list[ pt ].MOM[ mu ] = 0 ;
    }
    // equivalent Fourier modes in the -Pi -> Pi BZ
    list[ pt ].MOM[ ND-1 ] = pt > ( lt_2 - 1 ) ? -L0 + pt : pt ;
    list[ pt ].idx = pt ;

    // precompute all the twiddles along the t-direction
    double complex twiddles[ L0 ] ;
    int t ;
    double cache = 0.0 ;
    for( t = 0 ; t < L0 ; t++ ) {
      cache = ( ( pt ) * t * Latt.twiddles[ ND-1 ] ) ;
      twiddles[ t ] = cos( cache ) - I * sin( cache ) ;
    }

    // loop all polarisations
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	register double complex sumAA = 0.0 ;
	register double complex sumVV = 0.0 ;
	for( t = 0 ; t < L0 ; t++ ) {
	  sumAA += twiddles[ t ] * ctAA[ mu ][ nu ].mom[ 0 ].C[ t ] ;
	  sumVV += twiddles[ t ] * ctVV[ mu ][ nu ].mom[ 0 ].C[ t ] ;
	}
	cpAA[ pt ].PI[ mu ][ nu ] = sumAA ;
	cpVV[ pt ].PI[ mu ][ nu ] = sumVV ;
      }
    }
    //
  }
  // and return the veclist
  return list ;
}

// sum over every timeslice
static void
spatialzero_DFT( struct mcorr **corr ,
		 const struct PIdata *data )
{
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
	corr[ mu ][ nu ].mom[ 0 ].C[ t ] = sum ;
      }
    }
  }
  return ;
}

// time-moments approach
void
tmoments( const struct PIdata *AA ,
	  const struct PIdata *VV ,
	  const char *outfile ,
	  const current_type current ) 
{
  // set up the veclist
  int *tNMOM = malloc( sizeof( int ) ) ;
  const struct veclist *tlist = zero_veclist( tNMOM , ND-1 , GLU_FALSE ) ;

  // storage for the temporal data
  struct mcorr **ctAA = allocate_momcorrs( ND , ND , tNMOM[0] ) ;
  struct mcorr **ctVV = allocate_momcorrs( ND , ND , tNMOM[0] ) ;

  // sum over every timeslice
  spatialzero_DFT( ctAA , AA ) ;
  spatialzero_DFT( ctVV , VV ) ;

  // change our output files
  char strAA[ 256 ] , strVV[ 256 ] ;
  switch( current ) {
  case CONSERVED_LOCAL :
    sprintf( strAA , "%s.CALA.tcorr.bin" , outfile ) ;
    sprintf( strVV , "%s.CVLV.tcorr.bin" , outfile ) ;
    break ;
  case LOCAL_LOCAL :
    sprintf( strAA , "%s.LALA.tcorr.bin" , outfile ) ;
    sprintf( strVV , "%s.LVLV.tcorr.bin" , outfile ) ;
    break ;
  }

  // write out the t-correlators
  write_momcorr( strAA , (const struct mcorr**)ctAA , 
		 tlist , ND , ND , tNMOM ) ;
  write_momcorr( strVV , (const struct mcorr**)ctVV ,
		 tlist , ND , ND , tNMOM ) ;

  // storage for the momentum-space data
  struct PIdata *cpAA = malloc( L0 * sizeof( struct PIdata ) ) ;
  struct PIdata *cpVV = malloc( L0 * sizeof( struct PIdata ) ) ;

  // we now have a correlator C(t) which we Fourier transform
  // in the t-direction
  const struct veclist *list = DFT( cpAA , cpVV , 
				    (const struct mcorr**)ctAA , 
				    (const struct mcorr**)ctVV ) ;

  const int NMOM[ 1 ] = { L0 } ;

  // free the temporal correlators now
  free_momcorrs( ctAA , ND , ND , tNMOM[0] ) ;
  free_momcorrs( ctVV , ND , ND , tNMOM[0] ) ;

  // free the temporal ones
  free( tNMOM ) ; free( (void*)tlist ) ;

  // allocate momenta
  double *psq = malloc( NMOM[0] * sizeof( double ) ) ;
  double **p = malloc( NMOM[0] * sizeof( double* ) ) ;

  int i ;
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    p[i] = (double*)malloc( NMOM[0] * sizeof( double ) ) ;
  }

  // precompute momenta
  compute_p_psq( p , psq , list , NMOM[0] ) ;

  // change the output
  char str[ 256 ] ;
  sprintf( str , "%s.ptonly" , outfile ) ;

  momspace_data( cpAA , (const double **)p , psq , list , 
		 NMOM , str , current , AXIAL ) ;

  momspace_data( cpVV , (const double **)p , psq , list , 
		 NMOM , str , current , VECTOR ) ;

  // free the momentum list
  free( (void*)list ) ;

  // free psq and p
  free( psq ) ;
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    free( p[i] ) ;
  }
  free( p ) ;

  // free the DFTd data
  free( cpAA ) ;
  free( cpVV ) ;

  return ;
}
