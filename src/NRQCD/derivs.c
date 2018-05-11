/**
   @file derivs.c
   @brief derivative terms
 */
#include "common.h"

#include "matrix_ops.h"
#include "halfspinor_ops.h"

// computes G(x) - U^\dagger(x-mu) G(x-\mu)
static void
gradback( struct halfspinor *der ,
	  const struct halfspinor *S ,
	  const size_t t ,
	  const size_t mu )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    const size_t bck  = lat[i].back[mu] ;
    const size_t bck2 = lat[i + t*LCU].back[mu] ;
    size_t d , j ;
    for( d = 0 ; d < NS ; d++ ) {
      multabdag( (void*)der[i].D[d] , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
      for( j = 0 ; j < NCNC ; j++ ) {
	der[i].D[d][j] = S[i].D[d][j] - der[i].D[d][j] ;
      }
    }
  }
  return ;
}

// computes U G(x+\mu) - G(x)
static void
gradforw( struct halfspinor *der ,
	  const struct halfspinor *S ,
	  const size_t t ,
	  const size_t mu )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    const size_t idx = i + t*LCU ;
    const size_t fwd  = lat[ i ].neighbor[mu] ;
    size_t d , j ;
    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)der[i].D[d] , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
      // might write a vectorisation for this
      for( j = 0 ; j < NCNC ; j++ ) {
	der[i].D[d][j] -= S[i].D[d][j] ;
      }
    }
  }
  return ;
}

// same as grad, sum over all mu
void
grad_sq( struct halfspinor *der2 ,
	 const struct halfspinor *S ,
	 const size_t t )
{      
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    
    zero_colormatrix( der2[i].D[0] ) ;
    zero_colormatrix( der2[i].D[1] ) ;
    zero_colormatrix( der2[i].D[2] ) ;
    zero_colormatrix( der2[i].D[3] ) ;

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t idx = i + t*LCU ;

    size_t mu , d ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t fwd  = lat[ i ].neighbor[mu] ;
      const size_t bck  = lat[ i ].back[mu] ;
      const size_t bck2 = lat[idx].back[mu] ;
      
      for( d = 0 ; d < NS ; d++ ) {
	multab( (void*)A , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
	multabdag( (void*)B , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
	// A = A + B 
	add_mat( (void*)A , (void*)B ) ;
	// A = A - 2 F -> S[i].D[d]
	colormatrix_Saxpy( A , S[i].D[d], -2. ) ;
	// der2 = A + der2
	add_mat( (void*)der2[i].D[d] , (void*)A ) ;
      }
    }
  }
  return ;
}


// computes grad^2 of Stemp and puts it in der
void
grad2( struct halfspinor *der2 ,
       const struct halfspinor *S ,
       const size_t t ,
       const size_t mu )
{      
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t idx = i + t*LCU ;

    const size_t fwd  = lat[ i ].neighbor[mu] ;
    const size_t bck  = lat[ i ].back[mu] ;
    const size_t bck2 = lat[idx].back[mu] ;

    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)der2[i].D[d] , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
      multabdag( (void*)A , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
      // DER += B 
      add_mat( (void*)der2[i].D[d] , (void*)A ) ;
      // DER = DER - 2 F -> S[i].D[d]
      colormatrix_Saxpy( der2[i].D[d] , S[i].D[d] , -2. ) ;
    }
  }
  return ;
}

// computes derivative of S and puts it in der
void
grad( struct halfspinor *der ,
      const struct halfspinor *S ,
      const size_t t ,
      const size_t mu )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t idx = i + t*LCU ;
    const size_t fwd  = lat[ i ].neighbor[mu] ;
    const size_t bck  = lat[ i ].back[mu] ;
    const size_t bck2 = lat[idx].back[mu] ;

    size_t d , j ;
    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)der[i].D[d] , (void*)lat[idx].O[mu] , (void*)S[fwd].D[d] ) ;
      multabdag( (void*)A , (void*)lat[bck2].O[mu] , (void*)S[bck].D[d] ) ;
      // might write a vectorisation for this
      for( j = 0 ; j < NCNC ; j++ ) {
	der[i].D[d][j] = ( der[i].D[d][j] - A[j] )/2. ;
      }
    }
  }
  return ;
}

// computes improved derivative of S and puts it in der
void
grad_imp( struct halfspinor *der ,
	  const struct halfspinor *S ,
	  const size_t t ,
	  const size_t mu )
{
  struct halfspinor t1[ LCU ] , t2[ LCU ] ;
  
  // compute standard derivative
  grad( der , S , t , mu ) ;  
  // backward derivative on spinor
  gradback( t1 , S , t , mu ) ;
  // central one
  grad( t2 , t1 , t , mu ) ;
  // forward derivate on the previous result
  gradforw( t1 , t2 , t , mu ) ;

  // computes der = der - t3/6
  halfspinor_Saxpy( der , t1 , -1./6. ) ;
	   
  return ;
}

void
grad_sq_imp( struct halfspinor *der ,
	     const struct halfspinor *S ,
	     const size_t t )
{
  struct halfspinor t1[ LCU ] , t2[ LCU ] ;
  
  grad_sq( der , S , t ) ;

  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( t1 , S , t , mu ) ;
    grad2( t2 , t1 , t , mu ) ;
    halfspinor_Saxpy( der , t2 , -1./12. ) ;
  }
  
  return ;
}
