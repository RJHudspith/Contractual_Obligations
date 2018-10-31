/**
   @file rhoeta_contract.c
   @brief compute only the connected parts of a 1-flavor rho-eta dimeson
 */
#include "common.h"

#include "gammas.h"       // gt_Gdag_gt()
#include "contractions.h" // simple_meson_contract()
#include "spinor_ops.h"   // spinmul_atomic_right()

double complex
rhoeta_contract( const struct spinor S ,
		 const struct spinor Sadj ,
		 const size_t GSRC ,
		 const size_t GSNK ,
		 const struct gamma *GAMMAS )
{
  // accumulated sum
  register double complex sum = 0.0 ;

  // precompute the tilded gammas
  const struct gamma gt_GSNK_gt = gt_Gdag_gt( GAMMAS[ GSNK ] ,
					      GAMMAS[ GAMMA_T ] ) ;
  const struct gamma gt_G5_gt = gt_Gdag_gt( GAMMAS[ GAMMA_5 ] ,
					    GAMMAS[ GAMMA_T ] ) ;

  // products of traces #1
  sum += \
    simple_meson_contract( gt_G5_gt , Sadj , GAMMAS[ GAMMA_5 ] , S ) *
    simple_meson_contract( gt_GSNK_gt , Sadj , GAMMAS[ GSRC ] , S ) ;

  // products of traces #2
  sum += \
    simple_meson_contract( gt_GSNK_gt , Sadj , GAMMAS[ GAMMA_5 ] , S ) *
    simple_meson_contract( gt_G5_gt , Sadj , GAMMAS[ GSRC ] , S ) ;
  
  // single trace #1
  struct spinor S1 = Sadj , Stmp = S , S2 = Sadj ;
  gamma_mul_lr( &Stmp , GAMMAS[ GSRC ] , gt_GSNK_gt ) ;
  spinmul_atomic_right( &S1 , Stmp ) ;
  Stmp = S ;
  gamma_mul_lr( &Stmp , GAMMAS[ GAMMA_5 ] , gt_G5_gt ) ;
  spinmul_atomic_right( &S2 , Stmp ) ;
  sum += bilinear_trace( S2 , S1 ) ;

  // single trace #2
  S1 = Sadj ; Stmp = S ; S2 = Sadj ;
  gamma_mul_lr( &Stmp , GAMMAS[ GSRC ] , gt_G5_gt ) ;
  spinmul_atomic_right( &S1 , Stmp ) ;
  Stmp = S ;
  gamma_mul_lr( &Stmp , GAMMAS[ GAMMA_5 ] , gt_GSNK_gt ) ;
  spinmul_atomic_right( &S2 , Stmp ) ;
  sum += bilinear_trace( S2 , S1 ) ;
  
  return sum ;
}
