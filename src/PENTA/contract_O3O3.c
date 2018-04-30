/**
   @file contract_O3O3.c
   @brief contract the second baryon-meson operator
 */
#include "common.h"

#include "contract_O2O2.h"

void
contract_O3O3( struct spinmatrix *P ,
	       double complex **F ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS ,
	       const uint8_t **loc )
{
  contract_O2O2( P , F , U , S , D , B , OP1 , OP2 , GAMMAS , loc ) ;
  return ;
}
