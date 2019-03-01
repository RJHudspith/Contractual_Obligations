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
static void
project_parity( double complex *pos ,
		double complex *neg ,
		const struct spinmatrix P ,
		const struct gamma GT ,
		const size_t op )
{
  const double complex t1 = spinmatrix_trace( P.D ) ;
  const double complex t2 = gammaspinmatrix_trace( GT , P.D ) ;
  switch( op%PENTA_NOPS ) {
    // diquark-y ones are positive parity ops
    case 0 : case 1 : case 2 : case 3 : case 6 :
    *neg = 0.5 * ( t1 + t2 ) ;
    *pos = 0.5 * ( t1 - t2 ) ;
    break ;
    // baryon-meson are negative parity ops and so we flip parity as the
    // backward-propagating state is the one we match to the diquarky ones
  case 4 : case 5 : case 7 : case 8 :
    *neg = 0.5 * ( t1 - t2 ) ;
    *pos = 0.5 * ( t1 + t2 ) ;
    break ;
  }
  return ;
}

// just does the (ud)(us)\bar{b} contraction
// I stick the positive parity in the first PENTA_NOPS indices and the
// negative parity ones in the PENTA_NOPS -> 2*PENTA_NOPS*PENTA_NBLOCK*PENTA_NBLOCK indices
// op1 has positive parity and op2 and op3 have negative so there is
// some swapping around necessary with the results
// operator matrix is interesting: call DG1 the diquark with G1 gamma matrix
// and M1G* the first baryon-meson and M2G* the second for PENTA_NBLOCK = 2
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

  
  size_t idx ;
  for( idx = 0 ; idx < (PENTA_NBLOCK*PENTA_NBLOCK*PENTA_NOPS) ; idx++ ) {
    
    // indexing                  
    const size_t b1 = (idx/(PENTA_NOPS*PENTA_NBLOCK)) ;
    const size_t b2 = (idx/PENTA_NOPS)%PENTA_NBLOCK ;
    const size_t i  = (idx)%PENTA_NOPS ;

    // turn of the diquarky ones until I fix them if I ever do
    switch( i ) {
    case 0 : case 1 : case 2 : case 3 : case 6 :
      continue ;
    default :
      break ;
    }
    
    // do the contractions
    struct spinmatrix P ;
    contract[i]( &P , F , L , L , S , bwdH ,
		 GBLOCK[ b1 ] ,
		 GBLOCK[ b2 ] ,
		 GAMMAS ,
		 loc ) ;

    // get the map correct
    const size_t icol = i%3 , irow = i/3 ;
    const size_t idx = icol*PENTA_NBLOCK + irow*3*PENTA_NBLOCK*PENTA_NBLOCK + b2 + b1*PENTA_NBLOCK*3 ;
    
    project_parity( result + idx ,
		    result + idx + PENTA_NOPS*PENTA_NBLOCK*PENTA_NBLOCK ,
		    P , GAMMAS[ GAMMA_T ] , i ) ;
  }

  return SUCCESS ;
}
