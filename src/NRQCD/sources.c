/**
   @file sources.c
   @brief NRQCD sources
*/
#include "common.h"

#include "geometry.h"
#include "halfspinor_ops.h" // zero_halfspinor

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
void
initialise_source( struct halfspinor *S ,
		   const sourcetype source ,
		   const double twists[ ND ] ,
		   const size_t origin[ ND ] )
{
  size_t i ;
  int or[ NS ] = { 0 , 0 , 0 , 0 } ;
  for( i = 0 ; i < NS-1 ; i++ ) {
    or[i] = (int)origin[i] ;
  }

  // point source position
  const size_t idx = gen_site( or ) ;
  
  // initialise source position
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    // zero all the spinors
    zero_halfspinor( &S[i] ) ;

    switch( source ) {
    case POINT :
      if( i == idx ) {
	set_prop_to_constant( &S[ i ] , 1.0 ) ;
      }
      break ;
    case WALL :
      set_prop_to_constant( &S[ i ] , get_eipx( twists , i , ND-1 ) ) ;
      break ;
    }
  }
  return ;
}
