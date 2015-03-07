/**
   @file tmoments_PImunu.c
   @brief computation of \f$PI(q^2)\f$ where \f$q=(0,0,0,q_t)\f$

   TODO :: readers & writers for this and momspace PImunu data
 */

#include "common.h"

// DFT in time direction, cp is the DFT in the time direction
static void
DFT( struct PIdata *cp ,
     const struct PIdata *ct )
{
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
  return ;
}

// do the projection to get \Pi(q^2)
static void
projection( double *trans ,
	    double *longitudinal ,
	    const struct PIdata *cp )
{
  // projection
  const int lt_2 = L0 >> 1 ;
  int pt ;
#pragma omp parallel for private(pt)
  for( pt = -lt_2 ; pt < lt_2 ; pt++ ) {
    const int pidx = pt + lt_2 ;
    const double psq = 4.0 * sin( 0.5 * pt * Latt.twiddles[ ND - 1 ] ) * sin( 0.5 * pt * Latt.twiddles[ ND - 1 ] ) ;
    const double invpsq = ( psq == 0.0 ) ? 1.0 : 1.0 / psq ;
    double complex sumtrans = 0.0 , sumlong = 0.0 ;
    int mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      sumtrans += cp[ pidx ].PI[ mu ][ mu ] ;
      sumlong  += -cp[ pidx ].PI[ ND-1 ][ ND-1 ] ;
    }
    trans[ pidx ] = creal( sumtrans ) * invpsq / (double)( ND - 1 ) ;
    longitudinal[ pidx ] = creal( sumlong ) * invpsq / (double)( ND - 1 ) ;
  }
  return ;
}

void
tmoments( const struct PIdata *data ) 
{
  // storage for the temporal data
  struct PIdata *ct = malloc( L0 * sizeof( struct PIdata ) ) ;

  // storage for the momentum-space data
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

  DFT( cp , ct ) ;

  double *trans = malloc( L0 * sizeof( double ) ) ;
  double *longitudinal = malloc( L0 * sizeof( double ) ) ;

  projection( trans , longitudinal , cp ) ;

  free( trans ) ;
  free( longitudinal ) ;

  free( ct ) ;
  free( cp ) ;

  return ;
}
