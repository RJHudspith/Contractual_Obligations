/**
   @file grad.c
   @brief single derivatives
 */
#include "common.h"

#include "halfspinor_ops.h"  // Fmunu halfspinor
#include "matrix_ops.h"      // add_mat
#include "mmul.h"            // multab

// computes improved derivative of S and puts it in der
void
grad_imp_LCU( struct halfspinor *der ,
	      const struct halfspinor *S ,
	      const struct field *Fmunu ,
	      const size_t t ,
	      const size_t mu )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // some temporaries we need
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    
    const size_t Uidx = i + t*LCU ;    
    const size_t Ubck = lat[ Uidx ].back[mu] ;
    
    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Sbck2 = lat[ Sbck ].back[mu] ;
    
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
      // computes der = -( der - A )/2
      colormatrix_Sa_xmy( der[i].D[d] , A , 2./3. ) ;

      // compute U(x)U(x+\mu)S(x+2\mu)
      multab( (void*)A ,
	      (void*)Fmunu[i].O[6+2*mu] ,
	      (void*)S[ Sfwd2 ].D[d] ) ;
      // compute U(x-mu)U(x-2mu)S(x-2\mu)
      multabdag( (void*)B ,
		 (void*)Fmunu[i].O[7+2*mu],
		 (void*)S[ Sbck2 ].D[d] ) ;
      colormatrix_Sa_xmy( A , B , -1./12. ) ;

      add_mat( (void*)der[i].D[d] , (void*)A ) ;
    }
  }
  return ;
}

// computes improved derivative of Fmunu.S and puts it in der
void
grad_imp_FMUNU( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  // some temporaries we need
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    
  const size_t Uidx = i + t*LCU ;    
  const size_t Ubck = lat[ Uidx ].back[mu] ;
    
  const size_t Sfwd = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
  const size_t Sbck = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;
  
  // precompute the prop left-multiplied by the clover 
  struct halfspinor fwd , bwd , fwd2 , bwd2 ;
  Fmunu_halfspinor( &fwd , Fmunu[ Sfwd ].O[ Fmunu_idx ] , S[ Sfwd ] ) ;
  Fmunu_halfspinor( &bwd , Fmunu[ Sbck ].O[ Fmunu_idx ] , S[ Sbck ] ) ;
  Fmunu_halfspinor( &fwd2 , Fmunu[ Sfwd2 ].O[ Fmunu_idx ] , S[ Sfwd2 ] ) ;
  Fmunu_halfspinor( &bwd2 , Fmunu[ Sbck2 ].O[ Fmunu_idx ] , S[ Sbck2 ] ) ;
    
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    // computes der = U(x) S(x+\mu)
    multab( (void*)der -> D[d] ,
	    (void*)lat[ Uidx ].O[mu] ,
	    (void*)fwd.D[d] ) ;
    // computes A = U^\dagger (x-\mu ) S(x-\mu)
    multabdag( (void*)A ,
	       (void*)lat[ Ubck ].O[mu] ,
	       (void*)bwd.D[d] ) ;
    // computes der = -( der - A )/2
    colormatrix_Sa_xmy( der -> D[d] , A , 2./3. ) ;

    // compute U(x)U(x+\mu)S(x+2\mu)
    multab( (void*)A ,
	    (void*)Fmunu[i].O[6+2*mu] ,
	    (void*)fwd2.D[d] ) ;
    // compute U(x-mu)U(x-2mu)S(x-2\mu)
    multabdag( (void*)B ,
	       (void*)Fmunu[i].O[7+2*mu],
	       (void*)bwd2.D[d] ) ;
    colormatrix_Sa_xmy( A , B , -1./12. ) ;

    add_mat( (void*)der -> D[d] , (void*)A ) ;
  }
  return ;
}

// computes the improved derivative of S and multiplies on the left by Fmunu
void
FMUNU_grad_imp( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  // some temporaries we need
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    
  const size_t Uidx = i + t*LCU ;    
  const size_t Ubck = lat[ Uidx ].back[mu] ;
    
  const size_t Sfwd = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
  const size_t Sbck = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;

  struct halfspinor res ;
    
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    // computes der = U(x) S(x+\mu)
    multab( (void*)res.D[d] ,
	    (void*)lat[ Uidx ].O[mu] ,
	    (void*)S[ Sfwd ].D[d] ) ;
    // computes A = U^\dagger (x-\mu ) S(x-\mu)
    multabdag( (void*)A ,
	       (void*)lat[ Ubck ].O[mu] ,
	       (void*)S[ Sbck ].D[d] ) ;
    // computes der = -( der - A )/2
    colormatrix_Sa_xmy( res.D[d] , A , 2./3. ) ;

    // compute U(x)U(x+\mu)S(x+2\mu)
    multab( (void*)A ,
	    (void*)Fmunu[i].O[6+2*mu] ,
	    (void*)S[ Sfwd2 ].D[d] ) ;
    // compute U(x-mu)U(x-2mu)S(x-2\mu)
    multabdag( (void*)B ,
	       (void*)Fmunu[i].O[7+2*mu],
	       (void*)S[ Sbck2 ].D[d] ) ;
    colormatrix_Sa_xmy( A , B , -1./12. ) ;

    add_mat( (void*)res.D[d] , (void*)A ) ;
  }
  // res is the improved gradient and we left multiply by the gauge field
  Fmunu_halfspinor( der , Fmunu[i].O[ Fmunu_idx ]  , res ) ;
  
  return ;
}
