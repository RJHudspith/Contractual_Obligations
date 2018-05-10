/**
   @file derivs.c
   @brief derivative terms
 */
#include "common.h"

#include "matrix_ops.h"

// same as grad, sum over all mu
void
grad_sq( struct halfspinor *der2 ,
	 const struct halfspinor *S ,
	 const size_t t )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
      
  size_t i , d , j ;
  for( i = 0 ; i < LCU ; i++ ) {

    const size_t idx = i + t*LCU ;

    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t bck  = lat[ i ].neighbor[mu] ;
      const size_t fwd  = lat[ i ].back[mu] ;
      const size_t bck2 = lat[idx].back[mu] ;
      
      for( d = 0 ; d < NS ; d++ ) {
	multab( (void*)A , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
	multabdag( (void*)B , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
	for( j = 0 ; j < NCNC ; j++ ) {
	  der2[i].D[d][j] += ( A[j] + B[j] - 2*S[i].D[d][j] ) ;
	}
      }
    }
  }
  return ;
}


// computes grad^2 of Stemp and puts it in der, atomically added
void
grad2( struct halfspinor *der2 ,
       const struct halfspinor *S ,
       const size_t t ,
       const size_t mu )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
      
  size_t i , d , j ;
  for( i = 0 ; i < LCU ; i++ ) {

    const size_t idx = i + t*LCU ;

    const size_t bck  = lat[ i ].neighbor[mu] ;
    const size_t fwd  = lat[ i ].back[mu] ;
    const size_t bck2 = lat[idx].back[mu] ;

    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)A , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
      multabdag( (void*)B , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
      for( j = 0 ; j < NCNC ; j++ ) {
	der2[i].D[d][j] += ( A[j] + B[j] - 2*S[i].D[d][j] ) ;
      }
    }
  }
  return ;
}

// computes derivative of S and puts it in der, atomically added
void
grad( struct halfspinor *der ,
      const struct halfspinor *S ,
      const size_t t ,
      const size_t mu )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

  size_t i , d , j ;
  for( i = 0 ; i < LCU ; i++ ) {

    const size_t idx = i + t*LCU ;

    const size_t bck  = lat[ i ].neighbor[mu] ;
    const size_t fwd  = lat[ i ].back[mu] ;
    const size_t bck2 = lat[idx].back[mu] ;

    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)A , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
      multabdag( (void*)B , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
      for( j = 0 ; j < NCNC ; j++ ) {
	der[i].D[d][j] += ( A[j] - B[j] ) ;
      }
    }
  }
  return ;
}
