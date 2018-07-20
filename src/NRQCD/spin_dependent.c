/**
   @file spin_dependent.c
   @brief spin dependent terms in NRQCD action C_3, C_4, C_7 , C_8, C_9 

   Pauli matrix definition

   \sigma_x = | 0 1 |
              | 1 0 |

   \sigma_y = |  0 -I |
              |  I  0 |

   \sigma_z = | 1  0 |
              | 0 -1 |
 **/
#include "common.h"

#include "grad.h"            // FMUNU_grad_imp() and grad_imp_FMUNU()
#include "grad_2.h"          // gradsq_imp()
#include "halfspinor_ops.h"  // halfspinor_Saxpy and alike
#include "mmul.h"            // multab

// puts result in F -> S3
// computes \grad x E, which is the determinant of
//
// |  x   y   z  |
// | dx  dy  dz  |
// | Fxt Fyt Fzt |
//
// dy Fzt - dz Fyt for the x - direction
// dz Fxt - dx Fzt for the y - direction
// dx Fyt - dy Fxt for the z - direction
//
// For the terms that are E x \grad we have the reverse
//
// |  x   y   z  |
// | Fxt Fyt Fzt |
// | dx  dy  dz  |
//
// Fyt dz - Fzt dy Fyt for the x - direction
// Fzt dx - Fxt dz Fzt for the y - direction
// Fxt dy - Fyt dx Fxt for the z - direction
//
// I then dot in whichever sigma matrix we are using in this round
static void
sigma_dot_grad_x_E( struct halfspinor *H ,
		    const struct field *Fmunu ,
		    const struct halfspinor *S ,
		    const double U_0 ,
		    const size_t i ,
		    const size_t t ,
		    const size_t mu1 ,
		    const size_t mu2 ,
		    const uint8_t sigma_map[ NS ] ,
		    const uint8_t imap[ NS ] )
{
  // minus index map, much like the gamma matrix technology
  const uint8_t mimap[ NS ] = { (imap[0]+2)%4 , (imap[1]+2)%4 ,
				(imap[2]+2)%4 , (imap[3]+2)%4 } ;

  struct halfspinor res ;
  grad_imp_FMUNU( &res , S , Fmunu , U_0 , i , t , mu1 , 3+mu2 ) ;
  halfspinor_sigma_Saxpy( H , res , sigma_map , imap ) ;
  grad_imp_FMUNU( &res , S , Fmunu , U_0 , i , t , mu2 , 3+mu1 ) ;
  halfspinor_sigma_Saxpy( H , res , sigma_map , mimap ) ;

  FMUNU_grad_imp( &res , S , Fmunu , U_0 , i , t , mu1 , 3+mu2 ) ;
  halfspinor_sigma_Saxpy( H , res , sigma_map , imap ) ;
  FMUNU_grad_imp( &res , S , Fmunu , U_0 , i , t , mu2 , 3+mu1 ) ;
  halfspinor_sigma_Saxpy( H , res , sigma_map , mimap ) ;
  
  return ;
}


// does \sigma.( \grad x E - E x \grad ) S
static void
sigma_gradxE( struct halfspinor *H , 
	      const struct field *Fmunu ,
	      const struct halfspinor *S ,
	      const double U_0 ,
	      const size_t i ,
	      const size_t t )
{
  zero_halfspinor( H ) ;
  
  // First is the x - direction
  const uint8_t sigma_x[ NS ] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imapx[NS] = { 0 , 0 , 0 , 0 } ;
  sigma_dot_grad_x_E( H , Fmunu , S , U_0 , i , t , 1 , 2 , sigma_x , imapx ) ;
  // Second is the y - direction
  const uint8_t sigma_y[NS] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imapy[NS] = { 3 , 3 , 1 , 1 } ;
  sigma_dot_grad_x_E( H , Fmunu , S , U_0 , i , t , 2 , 0 , sigma_y , imapy ) ;
  // third is the z-direction
  const uint8_t sigma_z[NS] = { 0 , 1 , 2 , 3} ;
  const uint8_t imapz[NS] = { 0 , 0 , 2 , 2 } ;
  sigma_dot_grad_x_E( H , Fmunu , S , U_0 , i , t , 0 , 1 , sigma_z , imapz ) ;
  
  return ;
}

// 
void
term_C3( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C3 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C3 / ( 2. * pow( 2*NRQCD.M_0 , 2 ) ) ;

  struct halfspinor res ;
  sigma_gradxE( &res , Fmunu , S , NRQCD.U0 , i , t ) ;

  halfspinor_Saxpy( H , res , fac ) ;
  
  return ;
}

// sigma.B*G
void
term_C4( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C4 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C4 / ( 2. * NRQCD.M_0 ) ;

  struct halfspinor res ;
  sigmaB_halfspinor( &res , Fmunu[i] , S[i] ) ;

  halfspinor_Saxpy( H , res , fac ) ;

  return ;
}

// this term is { grad^2 , \sigma.B } G
void
term_C7( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C7 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C7 / ( pow( 2*NRQCD.M_0 , 3 ) ) ;

  struct halfspinor res ;
  sigmaB_gradsq_imp( &res , S , Fmunu , i , t ) ;
  halfspinor_Saxpy( H , res , fac ) ;
  
  gradsq_imp_sigmaB( &res , S , Fmunu , i , t ) ;
  halfspinor_Saxpy( H , res , fac ) ;
  
  return ;
}

// this term is even worse than c3 -> use yet another temporary
void
term_C8( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C8 ) < NRQCD_TOL ) return ;
  
  const double fac = -3.0 * NRQCD.C8 / ( 4. * pow( 2*NRQCD.M_0 , 4 ) ) ;

  size_t i ;
  // does \grad^2 \sigma.\grad.E.G
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {  
    sigma_gradxE( &F -> S1[i] , F -> Fmunu , F -> S , NRQCD.U0 , i , t ) ;
  }
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    struct halfspinor res ;
    gradsq_imp( &res , F -> S1 , F -> Fmunu , i , t ) ;
    #ifdef LEGACY_NRQCD_COMPARE
    // (in)correction factor for our derivative
    halfspinor_Saxpy( &res , F -> S1[i] , ( 1. - 1/(NRQCD.U0*NRQCD.U0) )/2. ) ;
    #endif
    halfspinor_Saxpy( &F -> H[i] , res , fac ) ;
  }
  // does \sigma.\grad.E \grad^2 G
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    gradsq_imp( &F -> S1[i] , F -> S , F -> Fmunu , i ,  t ) ;
    #ifdef LEGACY_NRQCD_COMPARE
    // (in)correction factor for our derivative
    halfspinor_Saxpy( &F -> S1[i] , F -> S[i] ,
		      ( 1. - 1/(NRQCD.U0*NRQCD.U0) )/2. ) ;
    #endif
  }
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    struct halfspinor res ;
    sigma_gradxE( &res , F -> Fmunu , F -> S1 , NRQCD.U0 , i , t ) ;
    halfspinor_Saxpy( &F -> H[i] , res , fac ) ; 
  }

  return ;
}

// term is i\sigma.(ExE + BxB)
// in Randy's code the coefficients are split up, I won't do that for the moment
// A and B have the E fields and C and D have the B - fields
void
term_C9EB( struct halfspinor *H ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C9EB ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C9EB / ( pow( 2*NRQCD.M_0 , 3 ) ) ;

  struct halfspinor t1 , t2 ;
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex D[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  size_t j ;
  // x direction is F_13 F_23 - F_23 F_13
  zero_halfspinor( &t1 ) ;
  multab( (void*)A , (void*)Fmunu[i].O[4] , (void*)Fmunu[i].O[5] ) ;
  multab( (void*)B , (void*)Fmunu[i].O[5] , (void*)Fmunu[i].O[4] ) ;
  multab( (void*)C , (void*)Fmunu[i].O[1] , (void*)Fmunu[i].O[2] ) ;
  multab( (void*)D , (void*)Fmunu[i].O[2] , (void*)Fmunu[i].O[1] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    t1.D[1][j] += ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
    t1.D[2][j] += ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
  }

  // y direction is F_23 F_03 - F_03 F_23
  multab( (void*)A , (void*)Fmunu[i].O[5] , (void*)Fmunu[i].O[3] ) ;
  multab( (void*)B , (void*)Fmunu[i].O[3] , (void*)Fmunu[i].O[5] ) ;
  multab( (void*)C , (void*)Fmunu[i].O[2] , (void*)Fmunu[i].O[0] ) ;
  multab( (void*)D , (void*)Fmunu[i].O[0] , (void*)Fmunu[i].O[2] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    t1.D[1][j] += -I * ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
    t1.D[2][j] += +I * ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
  }
    
  // z direction is F_03 F_13 - F_13 F_03
  multab( (void*)A , (void*)Fmunu[i].O[3] , (void*)Fmunu[i].O[4] ) ;
  multab( (void*)B , (void*)Fmunu[i].O[4] , (void*)Fmunu[i].O[3] ) ;
  multab( (void*)C , (void*)Fmunu[i].O[0] , (void*)Fmunu[i].O[1] ) ;
  multab( (void*)D , (void*)Fmunu[i].O[1] , (void*)Fmunu[i].O[0] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    t1.D[0][j] +=  ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
    t1.D[3][j] += -( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
  }
  // halfspinor multiply put into F -> S2
  halfspinor_multiply( &t2 , t1 , S[i] ) ;
  halfspinor_iSaxpy( H , t2 , fac ) ;

  return ;
}
