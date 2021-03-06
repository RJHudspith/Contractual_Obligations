/**
   @file sources.c
   @brief NRQCD sources
*/
#include "common.h"

#include "geometry.h"       // get_eipx()
#include "halfspinor_ops.h" // zero_halfspinor
#include "par_rng.h"        //
#include "quark_smear.h"    // quark_source_smear()

// set propagator to IdentityxConstant
static void
set_prop_to_constant( struct halfspinor *S ,
		      const double complex C )
{
#if (NS==4) && (NC==3)
  S -> D[0][0] = C ;
  S -> D[0][4] = C ;
  S -> D[0][8] = C ;
  S -> D[3][0] = C ;
  S -> D[3][4] = C ;
  S -> D[3][8] = C ;
#elif (NS==4) && (NC==2)
  S -> D[0][0] = C ;
  S -> D[0][3] = C ;
  S -> D[3][0] = C ;
  S -> D[3][3] = C ;
#else
  size_t d , c ;
  for( d = 0 ; d < NS/2 ; d++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      S -> D[d*(NS/2+1)][c*(NC+1)] = C ;
    }
  }
#endif
  return ;
}

// for the moment point is at (0,0,0,t)
int
initialise_source( struct halfspinor *S ,
		   struct halfspinor *S1 ,
		   const struct propagator prop ) 
{
  size_t i ;
  int or[ NS ] , flag = SUCCESS ;
  for( i = 0 ; i < ND-1 ; i++ ) {
    or[i] = (int)prop.origin[i] ;
  }
  or[ ND-1 ] = 0 ;

  // function pointer for stochastic sources
  double complex (*noise)( const uint32_t thread ) ;
  noise = Z2xZ2 ;

  switch( prop.Source.type ) {
  case POINT : break ;
  case WALL : break ;
  case Z2_WALL :
    #pragma omp single
    {
      flag = initialise_par_rng( NULL ) ;
    }
    break ;
  case Z3_WALL :
    noise = Z3 ;
    #pragma omp single
    {
      flag = initialise_par_rng( NULL ) ;
    }
    break ;
  case U1_WALL :
    noise = U1 ;
    #pragma omp single
    {
      flag = initialise_par_rng( NULL ) ;
    }
    break ;
  }

  // point source position
  const size_t idx = gen_site( or ) ;
  const size_t Z2_sub = ( prop.Source.Z2_spacing > 1 ?\
			  prop.Source.Z2_spacing/2 : 1 ) ;
    
  // initialise source position
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // figure out if we hit a sparse z2 site or not
    int n[ ND ] ;
    get_mom_2piBZ( n , i , ND-1 ) ;
    size_t mu , sum = 0 ;
    GLU_bool sparse = GLU_TRUE ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      // treats the x,y,z part of the source as a cardinal shift
      const size_t shift =
	(size_t)(n[mu]-(int)prop.origin[mu]+(int)Latt.dims[mu])%Latt.dims[mu] ;
      sum += shift ;
      if( shift%(Z2_sub) != 0 ) { sparse = GLU_FALSE ; }
    }
    if( sum%prop.Source.Z2_spacing != 0 ) { sparse = GLU_FALSE ; }

    #ifdef VERBOSE
    if( sparse == GLU_TRUE ) {
      printf( "%d %d %d is in sparse \n" , n[0] , n[1] , n[2] ) ;
    } else {
      printf( "%d %d %d is not \n" , n[0] , n[1] , n[2] ) ;
    }
    #endif
    
    // zero all the spinors
    zero_halfspinor( &S[i] ) ;

    switch( prop.Source.type ) {
    case POINT :
      if( i == idx ) {
	set_prop_to_constant( &S[ i ] , get_eipx( prop.mom_source , i , ND-1 ) ) ;
      }
      break ;
    case WALL :
      if( compute_rsq( i , ND-1 ) < prop.Source.boxsize ) {
	set_prop_to_constant( &S[ i ] ,
			      get_eipx( prop.mom_source , i , ND-1 ) ) ;
      }
      break ;
    case Z2_WALL :
    case Z3_WALL :
    case U1_WALL :
      if( sparse == GLU_TRUE ) {
	set_prop_to_constant( &S[ i ] , 
			      get_eipx( prop.mom_source , i , ND-1 )
			      * noise( get_CORR_thread() ) ) ;
      }
      break ;
    }
  }
  
  // call the smearing -> totally works with momentum sources and twisting
  // in fact, twisting + quark smearing == momentum smearing
  // although noone seems to have pointed out that you can just smear a
  // momentum source too and get more or less the same result
  // Also, doesn't need to be a point source or anything
  if( prop.Source.smear == QUARK ) {
    source_smear( S , S1 , prop.origin[ND-1] , prop.Source ) ;
  }
  
  // clean up the RNG
  switch( prop.Source.type ) {
  case POINT : break ;
  case WALL : break ;
  case Z2_WALL :
  case Z3_WALL :
  case U1_WALL :
    #pragma omp single
    {
      free_par_rng( ) ;
    }
    break ;
  }
    
  return flag ;
}
