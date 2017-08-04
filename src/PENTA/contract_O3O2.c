/**
   @file contract_O2O3.c
   @brief perform the contraction of our 2 baryon-meson operators
 */
#include "common.h"

#include "contract_O2O3.h"      // contract_O2O3()

// contraction is just the same as the O2O3 with S and D interchanged
void
contract_O3O2( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  contract_O2O3( P , U , S , D , B , GAMMAS ) ;
  return ;
}
