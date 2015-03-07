/**
   @file WardIdentity.c
   @brief check the ward identity

   TODO :: we should probably have a config-space check also right?
 */

#include "common.h"

// check the ward identity
void
compute_WI( const struct PIdata *data ,
	    const double **p )
{
  // momenta and correction
  double sum = 0.0 , sum2 = 0.0 ;
  int i ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double complex loc_sum = 0.0 , loc_sum2 = 0.0 ;
    int mu , nu ;
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

// perform ward identity correction
void
correct_WI( struct PIdata *data ,
	    const correction_dir corr_dir )
{
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {

    // compute Fourier modes
    int k[ ND ] = {} ;
    TwoPI_mpipi_momconv( k , i , ND ) ;

    // precompute correction factors
    double complex epi[ ND ] ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      const double cache = k[ mu ] * Latt.twiddles[ mu ] ;
      epi[ mu ] = cos( cache * 0.5 ) - I * sin( cache * 0.5 ) ;
    }
    
    // loop directions
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	switch( corr_dir ) {
	case CORR_MU : data[i].PI[mu][nu] *= epi[ mu ] ; break ;
	case CORR_NU : data[i].PI[mu][nu] *= epi[ nu ] ; break ;
	case CORR_MUpNU : data[i].PI[mu][nu] *= epi[ mu ] * epi[ nu ] ; break ;
	case UNCORR : break ;
	}
      }
    }
    //
  }
  return ;
}

// computes momenta
void
compute_p_psq( double **p ,
	       double *psq )
{
  int i ;
#pragma omp parallel for private(i)
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
    }
  }

  return ;
}
