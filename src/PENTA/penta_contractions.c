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
// assumes OP1 is the gamma matrix for the forward propagating state
// and OP2 is the gamma for the backward propagating state
void (*contract[9])( struct spinmatrix *P ,
		     double complex **F ,
		     const struct spinor U ,
		     const struct spinor D ,
		     const struct spinor S ,
		     const struct spinor B ,
		     const struct gamma OP1 ,
		     const struct gamma OP2 ,
		     const struct gamma *GAMMAS ,
		     const uint8_t **loc ) = {
  contract_O1O1 , contract_O1O2 , contract_O1O3 ,
  contract_O2O1 , contract_O2O2 , contract_O2O3 ,
  contract_O3O1 , contract_O3O2 , contract_O3O3 } ;

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
// negative parity ones in the PENTA_NOPS -> 2*PENTA_NOPS*PENTA_NBLOCK*PENTA_NBLOCK indices
// op1 has positive parity and op2 and op3 have negative so there is
// some swapping around necessary with the results
// operator matrix is interesting: call DG1 the diquark with G1 gamma matrix
// and M1G* the first baryon-meson and M2G* the secon for PENTA_NBLOCK = 2
//
// |(DG1 DG1^*)  (DG1 DG2^*)  | (DG1  M1G1^*) (DG1  M1G2^*) | (DG1  M2G1^*) (DG1  M2G2^*)|
// |(DG2 DG1^*)  (DG2 DG2^*)  | (DG2  M1G1^*) (DG2  M1G2^*) | (DG2  M2G1^*) (DG2  M2G2^*)|
// |(M1G1 DG1^*) (M1G1 DG2^*) | (M1G1 M1G1^*) (M1G1 M1G2^*) | (M1G1 M2G1^*) (M1G1 M2G2^*)|
// |(M1G2 DG1^*) (M1G2 DG2^*) | (M1G2 M1G1^*) (M1G2 M1G2^*) | (M1G2 M2G1^*) (M1G2 M2G2^*)|
// |(M2G1 DG1^*) (M2G1 DG2^*) | (M2G1 M1G1^*) (M2G1 M1G2^*) | (M2G1 M2G1^*) (M2G1 M2G2^*)|
// |(M2G2 DG1^*) (M2G2 DG2^*) | (M2G2 M1G1^*) (M2G2 M1G2^*) | (M2G2 M2G1^*) (M2G2 M2G2^*)|
//
// so it has similar block matrix structure as the TETRA contractions
int
pentas( double complex *result ,
	double complex **F ,
	const struct spinor L , 
	const struct spinor S ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS ,
	const uint8_t **loc )
{
#if PENTA_NBLOCK > 4
  fprintf( stderr , "[PENTA] compiled PENTA_NBLOCK greater than we allow\n" ) ;
  return FAILURE ;
#endif
  struct gamma GBLOCK[ 4 ] = { GAMMAS[ GAMMA_5 ] ,
			       GAMMAS[ IDENTITY ] ,
			       GAMMAS[ AT ] ,
			       GAMMAS[ GAMMA_T ] } ;

  size_t i , b1 , b2 ;
  
  // traces
  for( b1 = 0 ; b1 < PENTA_NBLOCK ; b1++ ) {
    for( b2 = 0 ; b2 < PENTA_NBLOCK ; b2++ ) {
	
      for( i = 0 ; i < PENTA_NOPS ; i++ ) {
	
	struct spinmatrix P ;
	contract[i]( &P , F , L , L , S , bwdH ,
		     GBLOCK[ b1 ] ,
		     GBLOCK[ b2 ] ,
		     GAMMAS ,
		     loc ) ;

	// get the map correct
	const size_t irow = i/3 ;
	const size_t idx = i*2 + irow*PENTA_NBLOCK*3 + b2 + b1*PENTA_NBLOCK*3 ;
	
	project_parity( result + idx ,
			result + idx + PENTA_NOPS*PENTA_NBLOCK*PENTA_NBLOCK ,
			P , GAMMAS[ GAMMA_T ] ) ;
      }
    }
  }

  return SUCCESS ;
}
