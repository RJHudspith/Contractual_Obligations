/**
   @file grad_4.c
   @brief fourth order gradient terms
 */
#include "common.h"

#include "matrix_ops.h"     // add_mat()
#include "mmul.h"           // multab() etc
#include "halfspinor_ops.h" // zero_halfspinor() etc


// I make this
// U(x)U(x+\mu)S(x+\mu) - 4U(x)S(x+\mu) + 6S(x)
// - 4U^\dagger(x-mu)S(x-\mu) + U^\dagger(x-\mu)U^\dagger(x-2\mu)S(x-2\mu)
void
grad4( struct halfspinor *der ,
       const struct halfspinor *S ,
       const struct field *Fmunu ,
       const size_t i ,
       const size_t t ,
       const size_t mu )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  
  const size_t Uidx = i + t*LCU ;
  
  const size_t Sfwd  = lat[ i ].neighbor[mu] ;
  const size_t Sfwd2 = lat[ Sfwd ].neighbor[mu] ;
  
  const size_t Sbck  = lat[ i ].back[mu] ;
  const size_t Sbck2 = lat[ Sbck ].back[mu] ;
  
  const size_t Ubck = lat[ Uidx ].back[mu] ;
  
  size_t d ;
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
    constant_mul_gauge( der -> D[d] , -4. , B ) ;
    // DER = DER - 6 F -> S[i].D[d]
    colormatrix_Saxpy( der -> D[d] , S[i].D[d] , 6. ) ;
  
    // the extra terms two steps away!
    multab( (void*)A , (void*)Fmunu[i].O[6+2*mu] ,
	    (void*)S[ Sfwd2 ].D[d] ) ;
    
    multabdag( (void*)B , (void*)Fmunu[i].O[7+2*mu] ,
	       (void*)S[ Sbck2 ].D[d] ) ;
    add_mat( (void*)A , (void*)B ) ;
    add_mat( (void*)der -> D[d] , (void*)A ) ;
  }
  return ;
}

// extra special \grad^2 \grad^2 code. This thing is a beast and pretty
// slow but it uses no temporaries
void
grad_sqsq( struct halfspinor *der2 ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex D[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex E[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  
  zero_halfspinor( der2 ) ;
  
  const size_t Uidx = i + t*LCU ;

  size_t mu , nu , d ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    const size_t Sfwd = lat[ i ].neighbor[mu] ;
    const size_t Sbck = lat[ i ].back[mu] ;
    const size_t Ubck = lat[ Uidx ].back[mu] ;

    const size_t S_PmPm = lat[ lat[i].neighbor[mu] ].neighbor[mu] ;
    const size_t S_MmMm = lat[ lat[i].back[mu] ].back[mu] ;
    
    for( d = 0 ; d < NS ; d++ ) {
      // computes A = U(x) S(x+\mu)
      multab( (void*)A , (void*)lat[ Uidx ].O[mu] , (void*)S[ Sfwd ].D[d] ) ;
      // computes B = U^\dag(x-\mu) S(x-\mu)
      multabdag( (void*)B , (void*)lat[ Ubck ].O[mu] , (void*)S[ Sbck ].D[d] ) ;
      // A = A + B 
      add_mat( (void*)A , (void*)B ) ;
      
      // there turns out to be 6 of these factor altogether and an overall factor of -2
      colormatrix_Saxpy( der2 -> D[d] , A , -12. );

      multab( (void*)A , (void*)Fmunu[i].O[6+2*mu] , (void*)S[ S_PmPm ].D[d] ) ;
      add_mat( (void*)der2 -> D[d] , (void*)A ) ;
      
      multabdag( (void*)A , (void*)Fmunu[i].O[7+2*mu] , (void*)S[ S_MmMm ].D[d] ) ;
      add_mat( (void*)der2 -> D[d] , (void*)A ) ;
    }
  }

  // beginning of the mu-nu loop
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    for( nu = 0 ; nu < ND-1 ; nu++ ) {

      if( mu == nu ) continue ;
      
      const size_t S_PmPn = lat[ lat[i].neighbor[mu] ].neighbor[nu] ;
      const size_t S_PmMn = lat[ lat[i].neighbor[mu] ].back[nu] ;
      const size_t S_MmPn = lat[ lat[i].back[mu] ].neighbor[nu] ;
      const size_t S_MmMn = lat[ lat[i].back[mu] ].back[nu] ;

      const size_t U_Pm   = lat[ Uidx ].neighbor[mu] ;
      const size_t U_Mm   = lat[ Uidx ].back[mu] ;
      const size_t U_PmMn = lat[ lat[ Uidx ].neighbor[mu] ].back[nu] ;
      const size_t U_MmMn = lat[ lat[ Uidx ].back[mu] ].back[nu] ;
      
      // precompute these things
      multab_suNC( (void*)B , (void*)lat[ Uidx ].O[mu] , (void*)lat[ U_Pm ].O[nu] ) ;
      multabdag_suNC( (void*)C , (void*)lat[ U_Mm ].O[mu] , (void*)lat[ U_Mm ].O[nu] ) ;
      multab_dag_suNC( (void*)D , (void*)lat[ Uidx ].O[mu] , (void*)lat[ U_PmMn ].O[nu]  ) ;
      multab_dagdag_suNC( (void*)E , (void*)lat[ U_Mm ].O[mu] , (void*)lat[ U_MmMn ].O[nu] ) ;
      
      for( d = 0 ; d < NS ; d++ ) {
	multab( (void*)A , (void*)B , (void*)S[ S_PmPn ].D[d] ) ;
	add_mat( (void*)der2 -> D[d] , (void*)A ) ;
	
	multab( (void*)A , (void*)C , (void*)S[ S_MmPn ].D[d] ) ;
	add_mat( (void*)der2 -> D[d] , (void*)A ) ;
	  
	multab( (void*)A , (void*)D , (void*)S[ S_PmMn ].D[d] ) ;
	add_mat( (void*)der2 -> D[d] , (void*)A ) ;
	
	multab( (void*)A , (void*)E , (void*)S[ S_MmMn ].D[d] ) ;
	add_mat( (void*)der2 -> D[d] , (void*)A ) ;
      }
    }
  }
  
  // can add this right at the end for slightly fewer instruction requests
  // it is pretty funny that the answer to life the universe and everything
  // is in here
  halfspinor_Saxpy( der2 , S[i] , 42. ) ;

  return ;
}
