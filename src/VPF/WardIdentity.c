/**
   @file WardIdentity.c
   @brief check the ward identity
 */

#include "common.h"

// computes momenta
void
compute_p_psq( double **p ,
	       double *psq ,
	       const struct veclist *list ,
	       const int NMOM )
{
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {
    // compute psq and p[i]
    psq[i] = 0.0 ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      // lattice momentum
      p[i][mu] = 2.0 * sin( 0.5 * list[i].MOM[ mu ] * Latt.twiddles[ mu ] ) ;
      psq[i] += p[i][mu] * p[i][mu] ;
    }
  }

  return ;
}

// check the ward identity
void
compute_WI( const struct PIdata *data ,
	    const double **p ,
	    const struct veclist *list ,
	    const int NMOM )
{
  // momenta and correction
  double sum = 0.0 , sum2 = 0.0 ;
  int i ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum2)
  for( i = 0 ; i < NMOM ; i++ ) {

    // we just look at the momenta we are keeping
    const int list_idx = list[ i ].idx ;

    register double complex loc_sum = 0.0 , loc_sum2 = 0.0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	// compute ward identities
        loc_sum  += p[i][mu] * data[ list_idx ].PI[mu][nu] ;
	loc_sum2 += data[ list_idx ].PI[mu][nu] * p[i][nu] ;
      }
    }

    sum  = sum  + creal( loc_sum ) * creal( loc_sum ) + \
                  cimag( loc_sum ) * cimag( loc_sum ) ;
    sum2 = sum2 + creal( loc_sum2 ) * creal( loc_sum2 ) + \
                  cimag( loc_sum2 ) * cimag( loc_sum2 ) ;
    // end of i loop
  }

  const double NORM = 1.0 / (double)( NMOM * ND ) ;
  printf( "|| p_{mu} Pi_{mu,nu} || :: %e\n" , sum * NORM ) ;
  printf( "|| PI_{mu,nu} p_{nu} || :: %e\n\n" , sum2 * NORM ) ;
}

// perform ward identity correction
void
correct_WI( struct PIdata *data ,
	    const correction_dir corr_dir ,
	    const struct veclist *list ,
	    const int NMOM )
{
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {

    // set the mapping from the momentum list
    const int list_idx = list[ i ].idx ;

    // precompute correction factors
    double complex epi[ ND ] ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      const double cache = list[ i ].MOM[ mu ] * Latt.twiddles[ mu ] ;
      epi[ mu ] = cos( cache * 0.5 ) - I * sin( cache * 0.5 ) ;
    }
    
    // loop directions
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	switch( corr_dir ) {
	case CORR_MU : data[ list_idx ].PI[mu][nu] *= epi[ mu ] ; break ;
	case CORR_NU : data[ list_idx ].PI[mu][nu] *= epi[ nu ] ; break ;
	case CORR_MUpNU : data[ list_idx ].PI[mu][nu] *= epi[ mu ] * epi[ nu ] ; break ;
	case UNCORR : break ;
	}
      }
    }
    //
  }
  return ;
}

// computes \delta_\mu V_\mu V_\nu 
void
WI_configspace( const struct PIdata *data ,
		const struct site *lat )
{
  double sum = 0.0 ;
  int i ;
#pragma omp parallel for private(i) reduction(+:sum) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double complex der = 0.0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	der += data[ i ].PI[mu][nu] - data[ lat[i].back[mu] ].PI[mu][nu] ;
      }
    }
    sum = sum + cabs( der ) ;
  }
  printf( "\n[VPF] config-space violation %e \n\n" , sum / ( double)LVOLUME ) ;
  return ;
}
