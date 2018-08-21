/**
   @file sources.c
   @brief NRQCD sources
*/
#include "common.h"

#include "geometry.h"       // get_eipx()
#include "grad_2.h"         // gradsq()
#include "halfspinor_ops.h" // zero_halfspinor
#include "par_rng.h"        //

#define WUPP_FAC (0.3)

// performs N iterations of Wuppertal smearing
static void
wuppertal_smear( struct halfspinor *S ,
		 struct halfspinor *S1 ,
		 const size_t t ,
		 const size_t Nsmear )
{
  const double fac = WUPP_FAC / ( 1 + WUPP_FAC * 2*(ND-1) ) ;
  size_t n ;
  for( n = 0 ; n < Nsmear ; n++ ) {
    size_t i ;
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      struct halfspinor der ;
      S1[i] = S[i] ;
      gradsq( &der , S , i , t ) ;
      halfspinor_Saxpy( &S1[i] , der , fac ) ;
    }
    // copy back
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      S[i] = S1[i] ;
    }
  }
  return ;
}

// set propagator to IdentityxConstant
static void
set_prop_to_constant( struct halfspinor *S1 ,
		      const double complex C )
{
#if (NS==4) && (NC==3)
  S1 -> D[0][0] = C ;
  S1 -> D[0][4] = C ;
  S1 -> D[0][8] = C ;
  S1 -> D[3][0] = C ;
  S1 -> D[3][4] = C ;
  S1 -> D[3][8] = C ;
#else
  size_t d , c ;
  for( d = 0 ; d < NS/2 ; d++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      S1 -> D[d*(NS/2+1)][c*(NC+1)] = C ;
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

  switch( prop.source ) {
  case Z2_WALL :
    #pragma omp single
    {
      flag = initialise_par_rng( NULL ) ;
    }
    break ;
  default : break ;
  }

  // point source position
  const size_t idx = gen_site( or ) ;
  
  // initialise source position
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    // zero all the spinors
    zero_halfspinor( &S[i] ) ;

    switch( prop.source ) {
    case POINT :
      if( i == idx ) {
	set_prop_to_constant( &S[ i ] , 1.0 ) ;
      }
      break ;
    case WALL :
      set_prop_to_constant( &S[ i ] ,
			    get_eipx( prop.mom_source , i , ND-1 ) ) ;
      break ;
    case Z2_WALL :
      set_prop_to_constant( &S[ i ] , 
			    get_eipx( prop.mom_source , i , ND-1 )
			    * Z2xZ2( get_CORR_thread() )
			    ) ;
      break ;
    }
  }

  // call the smearing?
  if( prop.smear == QUARK ) {
    wuppertal_smear( S , S1 , prop.origin[ND-1] , prop.Nsmear ) ;
  }
  
  // clean up the RNG
  switch( prop.source ) {
  case Z2_WALL :
    #pragma omp single
    {
      free_par_rng( ) ;
    }
    break ;
  default : break ;
  }
    
  return flag ;
}
