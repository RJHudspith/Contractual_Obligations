/**
   @file contract_O3O3.c
   @brief contract the second baryon-meson operator
 */
#include "common.h"

#include "contract_O2O2.h"

void
contract_O3O3( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  contract_O2O2( P , U , S , D , B , GAMMAS ) ;
  return ;
}
