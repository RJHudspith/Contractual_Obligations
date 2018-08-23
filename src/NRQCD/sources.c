/**
   @file sources.c
   @brief NRQCD sources
*/
#include "common.h"

#include "geometry.h"       // get_eipx()
#include "grad_2.h"         // gradsq()
#include "halfspinor_ops.h" // zero_halfspinor
#include "par_rng.h"        //

// performs N iterations of some smearing
// works like
// S = exp( asmear\grad^2 ) S
// by approximating the exp as the iteration
// S = ( 1 + asmear grad^2/Nsmear )^Nsmear S
// This is very like the C_0 term of NRQCD
static void
quark_smear( struct halfspinor *S ,
	     struct halfspinor *S1 ,
	     const size_t t ,
	     const struct source_info Source )
{
  const double fac = Source.smalpha/Source.Nsmear ;
  size_t n ;
  for( n = 0 ; n < Source.Nsmear ; n++ ) {
    size_t i ;
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      struct halfspinor der ;
      S1[i] = S[i] ;
      gradsq( &der , S , i , t ) ;
      halfspinor_Saxpy( &S1[i] , der , fac ) ;
    }
    // set S = S1, S1 = S
    #pragma omp single
    {
      struct halfspinor *t = S ;
      S  = S1 ;
      S1 = t ;
    }
  }
  return ;
}

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

  switch( prop.Source.type ) {
  case POINT : break ;
  case WALL : break ;
  case Z2_WALL :
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

    // figure out if we are hit a sparse z2 site
    int n[ ND ] ;
    get_mom_2piBZ( n , i , ND-1 ) ;
    size_t mu , sum = 0 ;
    GLU_bool sparse_Z2 = GLU_TRUE ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t shift =
	(size_t)(n[mu]-(int)prop.origin[mu]+(int)Latt.dims[mu])%Latt.dims[mu] ;
      sum += shift ;
      if( shift%(Z2_sub) != 0 ) { sparse_Z2 = GLU_FALSE ; }
    }
    if( sum%prop.Source.Z2_spacing != 0 ) { sparse_Z2 = GLU_FALSE ; }

    #ifdef VERBOSE
    if( sparse_Z2 == GLU_TRUE ) {
      printf( "%d %d %d is in sparse Z2 \n" , n[0] , n[1] , n[2] ) ;
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
      set_prop_to_constant( &S[ i ] ,
			    get_eipx( prop.mom_source , i , ND-1 ) ) ;
      break ;
    case Z2_WALL :
      if( sparse_Z2 == GLU_TRUE ) {
	set_prop_to_constant( &S[ i ] , 
			      get_eipx( prop.mom_source , i , ND-1 )
			      * Z2xZ2( get_CORR_thread() )
			      ) ;
      }
      break ;
    }
  }
  
  // call the smearing?
  if( prop.Source.smear == QUARK ) {
    quark_smear( S , S1 , prop.origin[ND-1] , prop.Source ) ;
  }
  
  // clean up the RNG
  switch( prop.Source.type ) {
  case POINT : break ;
  case WALL : break ;
  case Z2_WALL :
    #pragma omp single
    {
      free_par_rng( ) ;
    }
    break ;
  }
    
  return flag ;
}
