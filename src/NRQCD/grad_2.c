/**
   @file grad_sq.c
   @brief routines for the second order derivatives
 */
#include "common.h"

#include "matrix_ops.h"      // add_mat() etc
#include "mmul.h"            // multab() etc
#include "halfspinor_ops.h"  // zero_halfspinor() etc

// same as grad, sum over all mu
void
gradsq( struct halfspinor *der2 ,
	const struct halfspinor *S ,
	const size_t i ,
	const size_t t )
{      
  zero_halfspinor( der2 ) ;

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
      add_mat( (void*)der2 -> D[d] , (void*)A ) ;
    }
  }
  return ;
}

// same as grad, sum over all mu
void
grad_sq_LCU( struct halfspinor *der2 ,
	     const struct halfspinor *S ,
	     const size_t t )
{
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    zero_halfspinor( &der2[i] ) ;

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

// improved \grad^2
void
gradsq_imp( struct halfspinor *der ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  
  const size_t Uidx = i + t*LCU ;

  zero_halfspinor( der ) ;

  size_t mu , d ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
  
    const size_t Sfwd  = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck  = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
    const size_t Ubck = lat[ Uidx ].back[mu] ;
    
    for( d = 0 ; d < NS ; d++ ) {
      
      // computes der = U S(x+\mu)
      multab( (void*)B ,
	      (void*)lat[ Uidx ].O[mu] ,
	      (void*)S[ Sfwd ].D[d] ) ;
      // computes A = U^\dag(x-\mu) S(x-\mu)
      multabdag( (void*)A ,
		 (void*)lat[ Ubck ].O[mu] ,
		 (void*)S[ Sbck ].D[d] ) ;
      // DER = -4( B + A )
      add_mat( (void*)B , (void*)A ) ;
      constant_mul_gauge( A , 4./3. , B ) ;
      add_mat( (void*)der -> D[d] , (void*)A ) ;
      // DER = DER - 6 F -> S[i].D[d]
      colormatrix_Saxpy( der -> D[d] , S[i].D[d] , -5./2. ) ;
      
      // the extra terms two steps away!
      multab( (void*)A , (void*)Fmunu[i].O[6+2*mu] ,
	      (void*)S[ Sfwd2 ].D[d] ) ;
      
      multabdag( (void*)B , (void*)Fmunu[i].O[7+2*mu] ,
		 (void*)S[ Sbck2 ].D[d] ) ;
      add_mat( (void*)A , (void*)B ) ;
      
      colormatrix_Saxpy( (void*)der -> D[d] , (void*)A , -1./12. ) ;
    }
  }
  return ;
}
