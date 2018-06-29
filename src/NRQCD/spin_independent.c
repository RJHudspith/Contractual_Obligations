/**
   @file spin_independent.c
   @brief collection of terms in the NRQCD hamiltonian
 */
#include "common.h"

#include "grad.h"           // FMUNU_grad_imp()
#include "grad_2.h"         // gradsq
#include "grad_4.h"         // grad_4() and grad_sqsq()
#include "halfspinor_ops.h" // halfspinor_iSaxpy
#include "matrix_ops.h"     // add_mat()
#include "mmul.h"           // multab

// This part computes -C_0 * ( grad^2 / 2*M_0 ) S(x,t) 
void
term_C0( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C0 ) < NRQCD_TOL ) return ;

  const double fac = -NRQCD.C0 / ( 2. * NRQCD.M_0 ) ;

  struct halfspinor res ;
  gradsq( &res , S , i , t ) ;
  halfspinor_Saxpy( H , res , fac ) ;
  
  return ;
}

// these are more complicated ~ grad^2 ( grad^2 G )
void
term_C1_C6( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C1 ) < NRQCD_TOL &&
      fabs( NRQCD.C6 ) < NRQCD_TOL ) return ;

  // C6 has a factor of the power of the hamiltonian in it
  const double fac = -NRQCD.C1 / pow( 2*NRQCD.M_0 , 3 ) 
    -NRQCD.C6 / ( 4. * NRQCD.N * pow( 2*NRQCD.M_0 , 2 ) ) ;

  struct halfspinor res ;
  grad_sqsq( &res , S , Fmunu , i , t ) ;
  halfspinor_Saxpy( H , res , fac ) ;
  
  return ;
}

// this term is i( \grad.E - E.\grad )
void
term_C2( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C2 ) < NRQCD_TOL ) return ;

  // important to note the + sign for this parameter!
  const double fac = NRQCD.C2 / ( 2. * pow( 2*NRQCD.M_0 , 2 ) ) ;
  
  // first term is i\grad.E.S second is iE.\grad.S 
  struct halfspinor res ;
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    // does i \grad FMUNU G
    grad_imp_FMUNU( &res , S , Fmunu , NRQCD.U0 , i , t , mu , ND-1+mu ) ;
    halfspinor_iSaxpy( H , res , +fac ) ;
    // does -i FMUNU \grad G
    FMUNU_grad_imp( &res , S , Fmunu , NRQCD.U0 , i , t , mu , ND-1+mu ) ;
    halfspinor_iSaxpy( H , res , -fac ) ; 
  }

  return ;
}
// this term is (\grad)^2
void
term_C5( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C5 ) < NRQCD_TOL ) return ;

  //  important to note the sign on this parameter is +
  const double fac = NRQCD.C5 / ( 24. * NRQCD.M_0 ) ;

  struct halfspinor res ;
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {

    grad4( &res , S , Fmunu , i , t , mu ) ;

    halfspinor_Saxpy( H , res , fac ) ;
  }

  return ;
}

// this term is E.E + B.B
void
term_C10EB( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C10EB ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C10EB / ( pow( 2 * NRQCD.M_0 , 3 ) ) ;

  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  size_t d ;
  multab( (void*)A , (void*)Fmunu[i].O[0] , (void*)Fmunu[i].O[0] ) ;
  for( d = 1 ; d < (ND-1)*(ND-2) ; d++ ) {
    multab( (void*)B , (void*)Fmunu[i].O[d] , (void*)Fmunu[i].O[d] ) ;
    add_mat( (void*)A , (void*)B ) ;
  }
  struct halfspinor res ;
  colormatrix_halfspinor( (void*)res.D , (const void*)A , (const void*)S[i].D ) ;

  halfspinor_Saxpy( H , res , fac ) ;
  
  return ;
}

void
term_C11( struct NRQCD_fields *F ,
	  const size_t t ,
	  const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C11 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C11 / ( 24. * pow( NRQCD.N , 2 ) * pow( 2*NRQCD.M_0 , 3 ) ) ;

  grad_sq_LCU( F -> S1 , F -> S , t ) ;
  grad_sq_LCU( F -> S2 , F -> S1 , t ) ;
  grad_sq_LCU( F -> S1 , F -> S2 , t ) ;

  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    halfspinor_Saxpy( &F -> H[i] , F -> S1[i] , fac ) ;
  }
  
  return ;
}
