/**
   @file diquark_contraction.c
   @brief perform a diquark contraction
 */
#include "common.h"

#include "contractions.h"
#include "gammas.h"
#include "matrix_ops.h"
#include "spinmatrix_ops.h"
#include "spinor_ops.h"

// diquark is epsilon_{eab}( \psi_a^T C_GSRC \chi_b )
// contraction reads
// 
// Tr_S( P_aa C_GSRC C_bb C_GSNK ) - Tr_S( P_ab C_GSRC C_ba C_GSNK )
// where C_GSNK is tilded
double complex
diquark( const struct spinor S1 , 
	 const struct spinor S2 ,
	 const struct gamma C_GSRC , 
	 const struct gamma C_GSNK ) 
{
  struct spinor ST = transpose_spinor( S1 ) ;

  double complex D1[ NSNS ] __attribute__((aligned(ALIGNMENT))) ;
  double complex D2[ NSNS ] __attribute__((aligned(ALIGNMENT))) ;

  // trace out the color indices into spinmatrices
  colortrace_spinor( D1 , &ST ) ;
  colortrace_spinor( D2 , &S2 ) ;

  // multiply on the right with the gammas
  spinmatrix_gamma( D1 , C_GSRC ) ;
  spinmatrix_gamma( D2 , C_GSNK ) ;

  register double complex sum = 0.0 ;

  // note that there is an implicit sign in simple_meson_contract
  sum +=
    trace_prod_spinmatrices( D1 , D2 )     
    +simple_meson_contract( C_GSNK , ST , C_GSRC , S2 )
    ;

  return -sum ;
}
