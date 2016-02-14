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
	       const size_t NMOM )
{
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {
    // compute psq and p[i]
    psq[i] = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      // lattice momentum
      p[i][mu] = 2.0 * sin( 0.5 * Latt.twiddles[ mu ] * list[i].MOM[ mu ] ) ;
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
	    const size_t NMOM )
{
  // momenta and correction
  double sum = 0.0 , sum2 = 0.0 ;
  size_t i ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum2)
  for( i = 0 ; i < NMOM ; i++ ) {

    // we just look at the momenta we are keeping
    const size_t list_idx = list[ i ].idx ;

    register double complex loc_sum = 0.0 , loc_sum2 = 0.0 ;
    size_t mu , nu ;
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
	    const size_t NMOM )
{
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {

    // set the mapping from the momentum list
    const size_t list_idx = list[ i ].idx ;

    // precompute correction factors
    double complex epi[ ND ] ;
    size_t mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      const double cache = 0.5 * list[ i ].MOM[ mu ] * Latt.twiddles[ mu ] ;
      epi[ mu ] = cos( cache ) - I * sin( cache ) ;
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

// computes V( x ) - V( x - \mu )
void
WI_configspace_bwd( const struct PIdata *data ,
		    const struct site *lat )
{
  double sum = 0.0 ;
  size_t i ;
#pragma omp parallel for private(i) reduction(+:sum) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double complex der = 0.0 ;
    size_t mu , nu ;
    for( nu = 0 ; nu < ND ; nu++ ) {
      for( mu = 0 ; mu < ND ; mu++ ) {
	der += data[ i ].PI[mu][nu] - data[ lat[i].back[mu] ].PI[mu][nu] ;
      }
      sum = sum + cabs( der ) ;
    }
  }
  printf( "\n[VPF] backward config-space violation %e \n\n" , 
	  sum / ( double)LVOLUME ) ;
  return ;
}

// computes V( x + \mu ) - V( x )
void
WI_configspace_fwd( const struct PIdata *data ,
		    const struct site *lat )
{
  double sum = 0.0 ;
  size_t i ;
#pragma omp parallel for private(i) reduction(+:sum) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double complex der = 0.0 ;
    size_t mu , nu ;
    for( nu = 0 ; nu < ND ; nu++ ) {
      for( mu = 0 ; mu < ND ; mu++ ) {
	der += data[ lat[ i ].neighbor[ mu ] ].PI[mu][nu] - data[ i ].PI[mu][nu] ;
      }
      sum = sum + cabs( der ) ;
    }
  }
  printf( "\n[VPF] forward config-space violation %e \n\n" , 
	  sum / ( double)LVOLUME ) ;
  return ;
}

// computes V( x + \mu ) - V( x - \mu )
void
WI_configspace_sym( const struct PIdata *data ,
		    const struct site *lat )
{
  double sum = 0.0 ;
  size_t i ;
#pragma omp parallel for private(i) reduction(+:sum) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu , nu ;
    for( nu = 0 ; nu < ND ; nu++ ) {
      register double complex der = 0.0 ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	const size_t xmmu = lat[i].back[mu] ;
	const size_t xpmu = lat[i].neighbor[mu] ;
	der += ( data[ xpmu ].PI[mu][nu] - data[ xmmu ].PI[mu][nu] ) ;
      }
      sum = sum + cabs( der ) ;
    }
  }
  printf( "\n[VPF] symmetric config-space violation %e \n\n" , 
	  sum / ( double)LVOLUME ) ;
  return ;
}
