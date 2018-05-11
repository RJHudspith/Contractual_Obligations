/**
   @file sources.c
   @brief NRQCD sources
*/
#include "common.h"

#include "geometry.h"

// set propagator to IdentityxConstant
static void
set_prop_to_constant( struct halfspinor *S1 ,
		      const double complex C )
{
#if NS==4 && NC==3
  S1 -> D[0][0] = C ;
  S1 -> D[0][4] = C ;
  S1 -> D[0][8] = C ;
  S1 -> D[3][0] = C ;
  S1 -> D[3][4] = C ;
  S1 -> D[3][8] = C ;
  return ;
#else
  size_t d , c ;
  for( d = 0 ; d < NS/2 ; d++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      S1 -> D[d*(NS/2+1)][c*(NC+1)] = C ;
    }
  }
#endif
}

// for the moment point is at (0,0,0,t)
void
initialise_source( struct halfspinor *S ,
		   const sourcetype source ,
		   const size_t origin[ ND ]  )
{
  size_t i ;
  int or[ ND ] ;
  for( i = 0 ; i < ND ; i++ ) {
    or[i] = (int)origin[i] ;
  }
  
  // initialise source position
  switch( source ) {
  case POINT :
    fprintf( stdout , "[NRQCD] computing a point source at i = %zu\n" ,
	     gen_site( or ) ) ;
    set_prop_to_constant( &S[ gen_site( or ) ] , 1.0 ) ;
    break ;
    // wall source
  case WALL :
    fprintf( stdout , "[NRQCD] computing a Wall source at t_0 %zu\n" ,
	     origin[ND-1] ) ;

    for( i = 0 ; i < LCU ; i++ ) {
      set_prop_to_constant( &S[ i + origin[ND-1]*LCU ] , 1.0 ) ;
    }
    break ;
  }
  return ;
}
