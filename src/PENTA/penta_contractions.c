/**
   @file penta_contractions.c
   @brief pentaaquark contractions
 */
#include "common.h"

// contraction codes for 3x3 GEVP
#include "contract_O1O1.h"
#include "contract_O1O2.h" 
#include "contract_O1O3.h"

#include "contract_O2O1.h"      
#include "contract_O2O2.h" 
#include "contract_O2O3.h"

#include "contract_O3O1.h"
#include "contract_O3O2.h"  
#include "contract_O3O3.h"  

#include "spinmatrix_ops.h" // spinmatrix_trace()

// contraction function calls
void (*contract[9])( struct spinmatrix *P ,
		     const struct spinor U ,
		     const struct spinor D ,
		     const struct spinor S ,
		     const struct spinor B ,
		     const struct gamma *GAMMAS ) = {
  contract_O1O1 , contract_O1O2 , contract_O1O3 ,
  contract_O2O1 , contract_O2O2 , contract_O2O3 ,
  contract_O3O1 , contract_O3O2 , contract_O3O3 } ;

// get our idx from individual colors
size_t
idx( const size_t b , const size_t bp ,
     const size_t c , const size_t cp ,
     const size_t g , const size_t gp ,
     const size_t h , const size_t hp )
{
  return b + NC * ( bp + NC * ( c + NC * ( cp + NC * ( g + NC * ( gp + NC * ( h + NC * hp ) ) ) ) ) ) ;
}

// get our idx from individual colors shortened for 6 indices
size_t
idx2( const size_t b , const size_t bp ,
      const size_t c , const size_t cp ,
      const size_t g , const size_t gp )
{
  return b + NC * ( bp + NC * ( c + NC * ( cp + NC * ( g + NC * gp ) ) ) ) ;
}

// project to a particular parity
void
project_parity( double complex *pos ,
		double complex *neg ,
		const struct spinmatrix P ,
		const struct gamma GT )
{
  const double complex t1 = spinmatrix_trace( P.D ) ;
  const double complex t2 = gammaspinmatrix_trace( GT , P.D ) ;
  *neg = 0.5 * ( t1 + t2 ) ; // L4 baryon project
  *pos = 0.5 * ( t1 - t2 ) ; // L5 baryon project
  return ;
}

// just does the (ud)(us)\bar{b} contraction
// I stick the positive parity in the first PENTA_NOPS indices and the
// negative parity ones in the PENTA_NOPS -> 2*PENTA_NOPS indices
int
pentas( double complex *result ,
	const struct spinor L , 
	const struct spinor S ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS )
{
  // traces
  size_t i ;
  for( i = 0 ; i < PENTA_NOPS ; i++ ) {
    struct spinmatrix P ;
    contract[i]( &P , L , L , S , bwdH , GAMMAS ) ;
    project_parity( result + i , result + i + PENTA_NOPS ,
		    P , GAMMAS[ GAMMA_T ] ) ;
  }

  return SUCCESS ;
}
