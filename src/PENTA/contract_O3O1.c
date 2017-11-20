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
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS )
{
  contract_O2O1( P , U , S , D , B , OP1 , OP2 , GAMMAS ) ;
  return ;
}
