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
  //#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    const size_t Uidx = i + t*LCU ;
    const size_t Sbck = lat[i].back[mu] ;
    const size_t Ubck = lat[Uidx].back[mu] ;
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      multabdag( (void*)der[ i ].D[d] ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      colormatrix_Sa_xmy( der[i].D[d] , S[i].D[d] , -1. ) ;
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
  //#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    const size_t Uidx = i + t*LCU ;
    const size_t Sfwd  = lat[ i ].neighbor[mu] ;
    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      multab( (void*)der[ i ].D[d] ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      colormatrix_Sa_xmy( der[i].D[d] , S[i].D[d] , 1. ) ;
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
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    
    zero_colormatrix( der2[i].D[0] ) ;
    zero_colormatrix( der2[i].D[1] ) ;
    zero_colormatrix( der2[i].D[2] ) ;
    zero_colormatrix( der2[i].D[3] ) ;

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t Uidx = i + t*LCU ;

    size_t mu , d ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t Sfwd = lat[ i ].neighbor[mu] ;
      const size_t Sbck = lat[ i ].back[mu] ;
      const size_t Ubck = lat[ Uidx ].back[mu] ;
      
      for( d = 0 ; d < NS ; d++ ) {
	// computes A = U(x) S(x+\mu)
	multab( (void*)A ,
		(void*)lat[ Uidx ].O[mu] ,
		(void*)S[ Sfwd ].D[d] ) ;
	// computes B = U^\dag(x-\mu) S(x-\mu)
	multabdag( (void*)B ,
		   (void*)lat[ Ubck ].O[mu] ,
		   (void*)S[ Sbck ].D[d] ) ;
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

// computes mu-dependent grad^2 of Stemp and puts it in der
void
grad2( struct halfspinor *der2 ,
       const struct halfspinor *S ,
       const size_t t ,
       const size_t mu )
{      
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t Uidx = i + t*LCU ;

    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      // computes der = U S(x+\mu)
      multab( (void*)der2[ i ].D[d] ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      // computes A = U^\dag(x-\mu) S(x-\mu)
      multabdag( (void*)A ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      // DER += A
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
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

    const size_t Uidx = i + t*LCU ;
    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    size_t d ;
    for( d = 0 ; d < NS ; d++ ) {
      // computes der = U(x) S(x+\mu)
      multab( (void*)der[ i ].D[d] ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      // computes A = U^\dagger (x-\mu ) S(x-\mu)
      multabdag( (void*)A ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      // computes der = ( der - A )/2
      colormatrix_Sa_xmy( der[i].D[d] , A , 0.5 ) ;
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

  // compute standard derivative
  grad( der , S , t , mu ) ;

  // improvement doesn't work when threaded at the moment!!
#if 0
    struct halfspinor t1[ LCU ] , t2[ LCU ] ;
    // backward derivative on spinor
    gradback( t1 , S , t , mu ) ;
    // central one
    grad( t2 , t1 , t , mu ) ;
    // forward derivate on the previous result
    gradforw( t1 , t2 , t , mu ) ;
  }
// computes der = der - t3/6
//halfspinor_Saxpy( der , t1 , -1./6. ) ;
#endif	   
  return ;
}

void
grad_sq_imp( struct halfspinor *der ,
	     const struct halfspinor *S ,
	     const size_t t )
{
  // computes the standard grad^2 term we all know and love
  grad_sq( der , S , t ) ;

  // correction term
#if 0
  struct halfspinor t1[ LCU ] , t2[ LCU ] ;
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( t1 , S , t , mu ) ;
    grad2( t2 , t1 , t , mu ) ;
    halfspinor_Saxpy( der , t2 , -1./12. ) ;
  }
#endif
  return ;
}
