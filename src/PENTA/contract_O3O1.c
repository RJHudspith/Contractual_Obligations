/**
   @file contract_O3O1.c
   @brief contract the other Baryon-Meson operator with the diquarks
 */
#include "common.h"

#include "contract_O2O1.h"

void
contract_O3O1( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  contract_O2O1( P , U , S , D , B , GAMMAS ) ;
  return ;
}
