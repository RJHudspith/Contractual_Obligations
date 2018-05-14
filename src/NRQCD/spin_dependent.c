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

#include "derivs.h"
#include "matrix_ops.h"

#include "halfspinor_ops.h"

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
		    struct halfspinor *S3 ,
		    struct halfspinor *S2 ,
		    struct halfspinor *S1 ,
		    const struct field *Fmunu ,
		    const struct halfspinor *S ,
		    const size_t t ,
		    const size_t mu1 ,
		    const size_t mu2 ,
		    const size_t sigma_map[ NS ] ,
		    const uint8_t imap[ NS ] )
{
  const uint8_t mimap[ NS ] = { (imap[0]+2)%4 , (imap[1]+2)%4 ,
				(imap[2]+2)%4 , (imap[3]+2)%4 } ;
  size_t i ;
  // first term is \grad x E
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_halfspinor( &S1[i] , Fmunu[i].O[3+mu2] , S[i] ) ;
  }
  grad_imp( S3 , S2 , S1 , t , mu1 ) ;
  halfspinor_sigma_Saxpy( H , S3 , sigma_map , imap ) ;
    
  // and then do the minus term
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_halfspinor( &S1[i] , Fmunu[i].O[3+mu1] , S[i] ) ;
  }
  grad_imp( S3 , S2 , S1 , t , mu2 ) ;
  halfspinor_sigma_Saxpy( H , S3 , sigma_map , mimap ) ;

  ////////////////////////////////////////////////////////////////////////
  // do E x \grad here also -> Signs are correct as \grad x E = -E x \grad 
  grad_imp( S2 , S1 , S , t , mu1 ) ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_halfspinor( &S3[i] , Fmunu[i].O[3+mu2] , S2[i] ) ;
  }
  halfspinor_sigma_Saxpy( H , S3 , sigma_map , imap ) ;

  // and the minus term
  grad_imp( S2 , S1 , S , t , mu2 ) ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_halfspinor( &S3[i] , Fmunu[i].O[3+mu1] , S2[i] ) ;
  }
  halfspinor_sigma_Saxpy( H , S3 , sigma_map , mimap ) ;

  return ;
}

// computes \sigma.B G
static void
sigmaB_G( struct halfspinor *H ,
	  struct halfspinor *S1 ,
	  const struct field Fmunu ,
	  const struct halfspinor S ,
	  const double fac )
{
  // initialise sigma.B into S1
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    S1 -> D[0][j] =  Fmunu.O[2][j] ;
    S1 -> D[1][j] =  Fmunu.O[0][j] - I * Fmunu.O[1][j] ;
    S1 -> D[2][j] =  Fmunu.O[0][j] + I * Fmunu.O[1][j] ;
    S1 -> D[3][j] = -Fmunu.O[2][j] ;
  }

  // multiply S on the left by this result
  struct halfspinor res ;
  halfspinor_multiply( &res , *S1 , S ) ;
  
  colormatrix_Saxpy( H -> D[0] , res.D[0] , fac ) ;
  colormatrix_Saxpy( H -> D[1] , res.D[1] , fac ) ;
  colormatrix_Saxpy( H -> D[2] , res.D[2] , fac ) ;
  colormatrix_Saxpy( H -> D[3] , res.D[3] , fac ) ;
       
  return ;
}

// result is in F -> S3
// does \sigma_i ( \grad x E - E x \grad ) S
static void
sigma_gradxE( struct halfspinor *H , 
	      struct halfspinor *S3 ,
	      struct halfspinor *S2 ,
	      struct halfspinor *S1 ,
	      const struct field *Fmunu ,
	      const struct halfspinor *S ,
	      const size_t t )
{
  zero_halfspinor( H ) ;

  // does \sigma_x ( \grad_1 F_23 - \grad_2 F_13 )
  const size_t sigma_x[NS] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imapx[NS] = { 0 , 0 , 0 , 0 } ;
  sigma_dot_grad_x_E( H , S3 , S2 , S1 , Fmunu , S ,
		      t , 1 , 2 , sigma_x , imapx ) ;
  
  // does \sigma_y ( \grad_2 F_03 - \grad_0 F_23 )
  const size_t sigma_y[NS] = { 2 , 3 , 0 , 1 } ;
  const uint8_t imapy[NS] = { 3 , 3 , 1 , 1 } ;
  sigma_dot_grad_x_E( H , S3 , S2 , S1 , Fmunu , S ,
		      t , 2 , 0 , sigma_y , imapy ) ;
  
  // does \sigma_z ( \grad_0 F_13 - \grad_1 F_03 )
  const size_t sigma_z[NS] = { 0 , 1 , 2 , 3} ;
  const uint8_t imapz[NS] = { 0 , 0 , 2 , 2 } ;
  sigma_dot_grad_x_E( H , S3 , S2 , S1 , Fmunu , S ,
		      t , 0 , 1 , sigma_z , imapz ) ;
  return ;
}

// 
void
term_C3( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C3 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C3 / ( 2. * pow( 2*NRQCD.M_0 , 2 ) ) ;

  sigma_gradxE( F -> S4 , F -> S3 , F -> S2 , F -> S1 ,
		F -> Fmunu , F -> S , t ) ;

  halfspinor_Saxpy( F -> H , F -> S4 , fac ) ;
  
  return ;
}

// sigma.B*G
void
term_C4( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C4 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C4 / ( 2. * NRQCD.M_0 ) ;

  size_t i ;
#pragma omp for private(i)
  for ( i = 0  ; i < LCU ; i++ ) {
    sigmaB_G( &F -> H[i] , &F -> S1[i] , F -> Fmunu[i] , F -> S[i] , fac ) ;
  }
  return ;
}

// this term is { grad^2 , \sigma.B } G
void
term_C7( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C7 ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C7 / ( pow( 2*NRQCD.M_0 , 3 ) ) ;

  size_t i ;
  grad_sq_imp( F -> S3 , F -> S2 , F -> S1 , F -> S , t ) ;
#pragma omp for private(i)
  for ( i = 0  ; i < LCU ; i++ ) {
    sigmaB_G( &F -> H[i] , &F -> S2[i] , F -> Fmunu[i] , F -> S3[i] , fac ) ;
  }
  
  // second term sigma.B.\grad^2(G)
  zero_halfspinor( F -> S1 ) ;
#pragma omp for private(i)
  for ( i = 0  ; i < LCU ; i++ ) {
    sigmaB_G( &F -> S1[i] , &F -> S2[i] , F -> Fmunu[i] , F -> S[i] , fac ) ;
  }
  grad_sq_imp( F -> S4 , F -> S3 , F -> S2 , F -> S1 , t ) ;

  halfspinor_Saxpy( F -> H , F -> S4 , fac ) ;
  
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

  // does \sigma.\grad.E.\grad^2
  sigma_gradxE( F -> S4 , F -> S3 , F -> S2 , F -> S1 ,
		F -> Fmunu , F -> S , t ) ;
  grad_sq( F -> S1 , F -> S4 , t ) ;
  halfspinor_Saxpy( F -> H , F -> S1 , fac ) ;

  // does \grad^2.\sigma.\grad.E
  grad_sq( F -> S1 , F -> S , t ) ;
  sigma_gradxE( F -> S4 , F -> S5 , F -> S3 , F -> S2 ,
		F -> Fmunu , F -> S1 , t ) ;
  halfspinor_Saxpy( F -> H , F -> S4 , fac ) ;
  
  return ;
}

// term is i\sigma.(ExE + BxB)
// in Randy's code the coefficients are split up, I won't do that for the moment
void
term_C9EB( struct NRQCD_fields *F ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C9EB ) < NRQCD_TOL ) return ;
  
  const double fac = -NRQCD.C9EB / ( pow( 2*NRQCD.M_0 , 3 ) ) ;
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex B[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex C[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    double complex D[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    size_t j ;
    // x direction is F_13 F_23 - F_23 F_13
    multab( (void*)A , (void*)F -> Fmunu[i].O[4] , (void*)F -> Fmunu[i].O[5] ) ;
    multab( (void*)B , (void*)F -> Fmunu[i].O[5] , (void*)F -> Fmunu[i].O[4] ) ;
    multab( (void*)C , (void*)F -> Fmunu[i].O[1] , (void*)F -> Fmunu[i].O[2] ) ;
    multab( (void*)D , (void*)F -> Fmunu[i].O[2] , (void*)F -> Fmunu[i].O[1] ) ;
    for( j = 0 ; j < NCNC ; j++ ) {
      F -> S1[i].D[1][j] = I * ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
      F -> S1[i].D[2][j] = F -> S1[i].D[1][j] ;
    }
    // y direction is F_23 F_03 - F_03 F_23
    multab( (void*)A , (void*)F -> Fmunu[i].O[5] , (void*)F -> Fmunu[i].O[3] ) ;
    multab( (void*)B , (void*)F -> Fmunu[i].O[3] , (void*)F -> Fmunu[i].O[5] ) ;
    multab( (void*)C , (void*)F -> Fmunu[i].O[2] , (void*)F -> Fmunu[i].O[0] ) ;
    multab( (void*)D , (void*)F -> Fmunu[i].O[0] , (void*)F -> Fmunu[i].O[2] ) ;
    for( j = 0 ; j < NCNC ; j++ ) {
      F -> S1[i].D[1][j] =  ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
      F -> S1[i].D[2][j] = -( F -> S1[i].D[1][j] ) ;
    }
    // z direction is F_03 F_13 - F_13 F_03
    multab( (void*)A , (void*)F -> Fmunu[i].O[3] , (void*)F -> Fmunu[i].O[4] ) ;
    multab( (void*)B , (void*)F -> Fmunu[i].O[4] , (void*)F -> Fmunu[i].O[3] ) ;
    multab( (void*)C , (void*)F -> Fmunu[i].O[0] , (void*)F -> Fmunu[i].O[1] ) ;
    multab( (void*)D , (void*)F -> Fmunu[i].O[1] , (void*)F -> Fmunu[i].O[0] ) ;
    for( j = 0 ; j < NCNC ; j++ ) {
      F -> S1[i].D[0][j] = I * ( ( A[j] + C[j] ) - ( B[j] + D[j] ) ) ;
      F -> S1[i].D[3][j] = F -> S1[i].D[0][j] ;
    }
    // halfspinor multiply put into F -> S2
    halfspinor_multiply( &F -> S2[i] , F -> S1[i] , F -> S[i] ) ;
  }
  halfspinor_Saxpy( F -> H , F -> S2 , fac ) ;
  return ;
}
