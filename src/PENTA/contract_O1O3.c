/**
   @file contract_O1O3.c
   @brief perform the diquarks contraction with the second Baryon-Meson op
 */
#include "common.h"

#include "contract_O1O2.h"

// this is just the same contractions with the O1O2
void
contract_O1O3( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma *GAMMAS )
{
  contract_O1O2( P , U , S , D , B , GAMMAS ) ;
  return ;
}
