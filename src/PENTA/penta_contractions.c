/**
   @file penta_contractions.c
   @brief pentaaquark contractions
 */
#include "common.h"

#include "contract_O1O1.h"  // contract_O1O1()
#include "contract_O1O2.h"  // contract_O1O2()
#include "contract_O2O1.h"  // contract_O2O1()      
#include "contract_O2O2.h"  // contract_O2O2()

#include "spinmatrix_ops.h" // spinmatrix_trace()

// get our idx from individual colors
size_t
idx( const size_t b , const size_t bp ,
     const size_t c , const size_t cp ,
     const size_t g , const size_t gp ,
     const size_t h , const size_t hp )
{
  return b + NC * ( bp + NC * ( c + NC * ( cp + NC * ( g + NC * ( gp + NC * ( h + NC * hp ) ) ) ) ) ) ;
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
  struct spinmatrix P ;
  
  // contract the first operator
  contract_O1O1( &P , L , S , bwdH , GAMMAS ) ;
  project_parity( result + 0 , result + PENTA_NOPS ,
		  P , GAMMAS[ GAMMA_T ] ) ;

  // contract the Baryon-Meson - Diquarks operator
  contract_O1O2( &P , L , S , bwdH , GAMMAS ) ;
  project_parity( result + 1 , result + 1 + PENTA_NOPS ,
		  P , GAMMAS[ GAMMA_T ] ) ;
  
  // contract the Diquarks - Baryon-Meson operator
  contract_O2O1( &P , L , S , bwdH , GAMMAS ) ;
  project_parity( result + 2 , result + 2 + PENTA_NOPS ,
		  P , GAMMAS[ GAMMA_T ] ) ;
  
  // contract the Baryon-Meson operator
  contract_O2O2( &P , L , S , bwdH , GAMMAS ) ;
  project_parity( result + 3 , result + 3 + PENTA_NOPS ,
		  P , GAMMAS[ GAMMA_T ] ) ;

  return SUCCESS ;
}
