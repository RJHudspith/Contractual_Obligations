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
	      const double U_0 ,
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
      #ifdef LEGACY_NRQCD_COMPARE
      colormatrix_Sa_xmy( der[i].D[d] , A , (7+1/(U_0*U_0))/12. ) ;
      #else
      colormatrix_Sa_xmy( der[i].D[d] , A , 2./3. ) ;
      #endif
      
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
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  // some temporaries we need
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex D[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex E[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex F[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;

  const size_t Uidx = i + t*LCU ;    
  const size_t Ubck = lat[ Uidx ].back[mu] ;
    
  const size_t Sfwd = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
  const size_t Sbck = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;

  // precompute some common terms here -> products of links and clovers
  multab( (void*)C , (void*)lat[ Uidx ].O[mu] ,
	  (void*)Fmunu[ Sfwd ].O[ Fmunu_idx ] ) ;
  multabdag( (void*)D , (void*)lat[ Ubck ].O[mu] ,
	     (void*)Fmunu[ Sbck ].O[ Fmunu_idx ] ) ;
  multab( (void*)E , (void*)Fmunu[i].O[6+2*mu] ,
	  (void*)Fmunu[ Sfwd2 ].O[ Fmunu_idx ] ) ;
  multabdag( (void*)F , (void*)Fmunu[i].O[7+2*mu] ,
	     (void*)Fmunu[ Sbck2 ].O[ Fmunu_idx ] ) ;

  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    // computes der = U(x) S(x+\mu)
    multab( (void*)der -> D[d] , (void*)C , (void*)S[ Sfwd ].D[d] ) ;
    // computes A = U^\dagger (x-\mu ) S(x-\mu)
    multab( (void*)A , (void*)D , (void*)S[ Sbck ].D[d] ) ;
    // computes der = ( der - A )/2
    #ifdef LEGACY_NRQCD_COMPARE
    colormatrix_Sa_xmy( der -> D[d] , A , ( 7 + 1/(U_0*U_0) )/12. ) ;
    #else
    colormatrix_Sa_xmy( der -> D[d] , A , 2/3. ) ;
    #endif
    
    // compute U(x)U(x+\mu)S(x+2\mu)
    multab( (void*)A , (void*)E , (void*)S[ Sfwd2 ].D[d] ) ;
    // compute U(x-mu)U(x-2mu)S(x-2\mu)
    multab( (void*)B , (void*)F , (void*)S[ Sbck2 ].D[d] ) ;
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
		const double U_0 ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx )
{
  // some temporaries we need
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex D[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  
  const size_t Uidx = i + t*LCU ;    
  const size_t Ubck = lat[ Uidx ].back[mu] ;
    
  const size_t Sfwd = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
    
  const size_t Sbck = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;

  struct halfspinor res ;

  // it turns out a little faster to do this than call multabdag
  dagger_gauge( (void*)C , (void*)lat[ Ubck ].O[mu] ) ;
  dagger_gauge( (void*)D , (void*)Fmunu[i].O[7+2*mu] ) ;

  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    // computes der = U(x) S(x+\mu)
    multab( (void*)res.D[d] ,
	    (void*)lat[ Uidx ].O[mu] ,
	    (void*)S[ Sfwd ].D[d] ) ;
    // computes A = U^\dagger (x-\mu ) S(x-\mu)
    multab( (void*)A , (void*)C , (void*)S[ Sbck ].D[d] ) ;
    // computes res = 2/3. * ( res - A )
    #ifdef LEGACY_NRQCD_COMPARE
    colormatrix_Sa_xmy( res.D[d] , A , ( 7 + 1/(U_0*U_0) )/12. ) ;
    #else
    colormatrix_Sa_xmy( res.D[d] , A , 2/3. ) ;
    #endif
 
    // compute U(x)U(x+\mu)S(x+2\mu)
    multab( (void*)A , (void*)Fmunu[i].O[6+2*mu] , (void*)S[ Sfwd2 ].D[d] ) ;
    // compute U(x-mu)U(x-2mu)S(x-2\mu)
    multab( (void*)B , (void*)D , (void*)S[ Sbck2 ].D[d] ) ;
    // computes A = -1/12. * ( A - B )
    colormatrix_Sa_xmy( A , B , -1./12. ) ;

    // adds A to res
    add_mat( (void*)res.D[d] , (void*)A ) ;
  }
  // res is the improved gradient and we left multiply by the gauge field
  Fmunu_halfspinor( der , Fmunu[i].O[ Fmunu_idx ]  , res ) ;
  
  return ;
}
