/**
   @file contract_O1O3.c
   @brief perform the diquarks contraction with the second Baryon-Meson op
 */
#include "common.h"

#include "contract_O1O2.h"

// this is just the same contractions with the O1O2
void
contract_O1O3( struct spinmatrix *P ,
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
  contract_O1O2( P , F , U , S , D , B , OP1 , OP2 , GAMMAS , loc ) ;
  return ;
}
