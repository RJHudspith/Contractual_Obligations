/**
   @file diquark_contraction.c
   @brief perform a diquark contraction
 */
#include "common.h"

#include "spinmatrix_ops.h"
#include "tetra_contractions.h"

// diquark is ( \psi^T C_GSRC \chi )
// contraction reads
// 
// Tr_S( P_ac C_GSRC C_bd C_GSNK )
// where C_GSNK is tilded
int
diquark( double complex *result , 
	 struct spinor S1 , 
	 struct spinor S2 ,
	 const struct gamma C_GSRC , 
	 const struct gamma C_GSNK ) 
{
  // flattened number of colors we loop over
  const size_t Nco = NCNC * NCNC ;

  // index loops
  size_t abcd ;

  // temporaries blocks of spinmatrix data
  struct block *C1 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco * sizeof( struct block ) ) != 0 ) {
    goto memfree ;
  }

  // precompute usual spinor  block
  precompute_block( C1 , S1 , C_GSRC , S2 , C_GSNK ) ;

  // loop all possible color combinations
  for( abcd = 0 ; abcd < Nco ; abcd++ ) {
    result[ abcd ] = spinmatrix_trace( C1[ abcd ].M ) ;
  }

  free( C1 ) ;

  return SUCCESS ;

 memfree :
  free( C1 ) ;
  return sqrt(-1) ;
}
