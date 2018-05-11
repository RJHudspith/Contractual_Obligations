/**
   @file spin_independent.c
   @brief collection of terms in the NRQCD hamiltonian
 */
#include "common.h"

#include "derivs.h"
#include "matrix_ops.h"
#include "halfspinor_ops.h"

// This part computes -C_0 * ( grad^2 / 2*M_0 ) S(x,t) 
void
term_C0( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C0 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C0 / ( 2. * NRQCD.M_0 ) ;

  grad_sq( F -> S1 , F -> S , t ) ;

  halfspinor_Saxpy( F -> H , F -> S1 , fac ) ;

  return ;
}

// these are more complicated ~ grad^2 ( grad^2 G )
void
term_C1_C6( struct NRQCD_fields *F ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C1 ) < 1E-12 && fabs( NRQCD.C6 ) < 1E-12) return ;

  // C6 has a factor of the power of the hamiltonian in it
  const double fac = -NRQCD.C1 / pow( 2*NRQCD.M_0 , 3 ) -
    NRQCD.C6 / ( 4. * NRQCD.N * pow( 2*NRQCD.M_0 , 2 ) ) ;

  grad_sq( F -> S1 , F -> S , t ) ;
  grad_sq( F -> S2 , F -> S1 , t ) ;
  
  halfspinor_Saxpy( F -> H , F -> S2 , fac ) ;
  
  return ;
}

// this term is i( \grad.E - E.\grad )
void
term_C2( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C2 ) < 1E-12 ) return ;
  
  const double fac = NRQCD.C2 / ( 2. * pow( 2*NRQCD.M_0 , 2 ) ) ;
  
  size_t mu , i ;
  // first term is \grad.E
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      colormatrix_halfspinor( &F -> S1[i] ,
			      F -> Fmunu[i].O[ (ND-1) + mu ] ,
			      F -> S[i] ) ;
    }
    // accumulate derivative into S2
    grad_imp( F -> S2 , F -> S1 , t , mu ) ;
    halfspinor_iSaxpy( F -> H , F -> S2 , fac ) ;

    // call derivative into S1
    grad_imp( F -> S3 , F -> S , t , mu ) ;
    // E * prop goes into S2
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      halfspinor_colormatrix( &F -> S4[i] ,
			      F -> S3[i] ,
			      F -> Fmunu[i].O[ (ND-1) + mu ] ) ;
    }
    // update the Hamiltonian
    halfspinor_iSaxpy( F -> H , F -> S3 , -fac ) ;
  }

  return ;
}

// this term is (\grad)^2
void
term_C5( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C5 ) < 1E-12 ) return ;
  
  const double fac = NRQCD.C5 / ( 24. * NRQCD.M_0 ) ;

  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S1 , F -> S , t , mu ) ;
    grad2( F -> S2 , F -> S1 , t , mu ) ;
    halfspinor_Saxpy( F -> H , F -> S2 , fac ) ;
  }
  
  return ;
}

// this term is E.E + B.B
void
term_C10EB( struct NRQCD_fields *F ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C10EB ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C10EB / ( pow( 2 * NRQCD.M_0 , 3 ) ) ;

  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    size_t d ;
    multab( (void*)A , (void*)F -> Fmunu[i].O[0] , (void*)F -> Fmunu[i].O[0] ) ;
    for( d = 1 ; d < (ND-1)*(ND-2) ; d++ ) {
      multab( (void*)B , (void*)F -> Fmunu[i].O[d] , (void*)F -> Fmunu[i].O[d] ) ;
      add_mat( (void*)A , (void*)B ) ;
    }
    colormatrix_halfspinor( &F -> S1[i] , A , F -> S[i] ) ;
  }
  halfspinor_Saxpy( F -> H , F -> S1 , fac ) ;
  
  return ;
}

void
term_C11( struct NRQCD_fields *F ,
	  const size_t t ,
	  const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C11 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C11 / ( 192. * pow( NRQCD.N , 2 ) * pow( NRQCD.M_0 , 3 ) ) ;
  
  grad_sq( F -> S1 , F -> S , t ) ;
  grad_sq( F -> S2 , F -> S1 , t ) ;
  grad_sq( F -> S3 , F -> S2 , t ) ;

  halfspinor_Saxpy( F -> H , F -> S3 , fac ) ;
}
