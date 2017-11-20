/**
   @file contract_O2O2.c
   @brief contract the Baryon-Meson operators
 */
#include "common.h"

#include "bar_contractions.h" // bar_contract_site()
#include "contractions.h"     // gamma_mul_r()
#include "gammas.h"           // CGmu()

// baryon-meson contraction
void
contract_O2O2( struct spinmatrix *P ,
	       const struct spinor U ,
	       const struct spinor D ,
	       const struct spinor S ,
	       const struct spinor B ,
	       const struct gamma OP1 ,
	       const struct gamma OP2 ,
	       const struct gamma *GAMMAS )
{
  // usual gamma stuff
  const struct gamma G5    = OP1 ;
  const struct gamma C1    = CGmu( OP1 , GAMMAS ) ;
  
  const struct gamma t2t   = gt_Gdag_gt( OP2 , GAMMAS[ GAMMA_T ] ) ;
  const struct gamma tC2t  = gt_Gdag_gt( CGmu( OP2 , GAMMAS ) ,
					 GAMMAS[ GAMMA_T ] ) ;

  // allocat the term matrix
  double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
  term[0] = malloc( NSNS * sizeof( double complex ) ) ;
  term[1] = malloc( NSNS * sizeof( double complex ) ) ;

  // zero the term
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NSNS ; d1++ ) {
    term[0][d1] = term[1][d1] = 0.0 ;
  }

  baryon_contract_site( term , U , U , D , C1 , tC2t ) ;
  
  const double complex T =
    -simple_meson_contract( t2t , B , G5 , S ) ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) { 
      P -> D[d1][d2] =
	( term[0][ d2 + NS * d1 ] + term[1][ d2 + NS * d1 ] ) * T ;
    }
  }
  
  // free the term temporary
  free( term[0] ) ;
  free( term[1] ) ;
  free( term ) ;

  return ;
}
