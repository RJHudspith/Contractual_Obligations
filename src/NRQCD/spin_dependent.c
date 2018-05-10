/**
   @file spin_dependent.c
   @brief spin dependent terms in NRQCD action C_3, C_4, C_7 , C_8, C_9 ..
 **/
#include "common.h"

#include "derivs.h"
#include "matrix_ops.h"

// set s-matrix to zero
static void
set_S_to_zero( struct halfspinor *S )
{
  size_t i , j ;
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      S[i].D[0][j] = 0.0 ;
      S[i].D[1][j] = 0.0 ;
      S[i].D[2][j] = 0.0 ;
      S[i].D[3][j] = 0.0 ;
    }
  }
  return ;
}

static void
update_H( struct halfspinor *H ,
	  const double fac ,
	  const struct halfspinor *S )
{
  size_t i , j ;
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      H[i].D[0][j] += fac * S[i].D[0][j] ;
      H[i].D[1][j] += fac * S[i].D[1][j] ;
      H[i].D[2][j] += fac * S[i].D[2][j] ;
      H[i].D[3][j] += fac * S[i].D[3][j] ;
    }
  }
}

// puts result in F -> S3
static void
sigma_dot_grad_x_E( struct halfspinor *S3 ,
		    struct halfspinor *S2 ,
		    struct halfspinor *S1 ,
		    struct field *Fmunu ,
		    const struct halfspinor *S ,
		    const size_t t ,
		    const size_t Fidx1 ,
		    const size_t Fidx2 ,
		    const size_t mu1 ,
		    const size_t mu2 ,
		    const size_t sigma_map[ NS ] ,
		    const double complex mult_map[ NS ] )
{
  size_t i , j ;
  for( i = 0 ; i < LCU ; i++ ) {
    multab( (void*)S1[i].D[0] , (void*)Fmunu[i].O[Fidx1] , (void*)S[i].D[0] ) ;
    multab( (void*)S1[i].D[1] , (void*)Fmunu[i].O[Fidx1] , (void*)S[i].D[1] ) ;
    multab( (void*)S1[i].D[2] , (void*)Fmunu[i].O[Fidx1] , (void*)S[i].D[2] ) ;
    multab( (void*)S1[i].D[3] , (void*)Fmunu[i].O[Fidx1] , (void*)S[i].D[3] ) ;
  }
  set_S_to_zero( S2 ) ;
  grad( S2 , S1 , t , mu1 ) ;
  // multiply by sigma
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      S3[i].D[0][j] += mult_map[0] * S2[i].D[ sigma_map[0] ][j] ;
      S3[i].D[1][j] += mult_map[1] * S2[i].D[ sigma_map[1] ][j] ;
      S3[i].D[2][j] += mult_map[2] * S2[i].D[ sigma_map[2] ][j] ;
      S3[i].D[3][j] += mult_map[3] * S2[i].D[ sigma_map[3] ][j] ;
    }
  }
  // do the minus term
  for( i = 0 ; i < LCU ; i++ ) {
    multab( (void*)S1[i].D[0] , (void*)Fmunu[i].O[Fidx2] , (void*)S[i].D[0] ) ;
    multab( (void*)S1[i].D[1] , (void*)Fmunu[i].O[Fidx2] , (void*)S[i].D[1] ) ;
    multab( (void*)S1[i].D[2] , (void*)Fmunu[i].O[Fidx2] , (void*)S[i].D[2] ) ;
    multab( (void*)S1[i].D[3] , (void*)Fmunu[i].O[Fidx2] , (void*)S[i].D[3] ) ;
  }
  set_S_to_zero( S2 ) ;
  grad( S2 , S1 , t , mu2 ) ;
  // multiply by -sigma
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      S3[i].D[0][j] -= mult_map[0] * S2[i].D[ sigma_map[0] ][j] ;
      S3[i].D[1][j] -= mult_map[1] * S2[i].D[ sigma_map[1] ][j] ;
      S3[i].D[2][j] -= mult_map[2] * S2[i].D[ sigma_map[2] ][j] ;
      S3[i].D[3][j] -= mult_map[3] * S2[i].D[ sigma_map[3] ][j] ;
    }
  }
  ////////////////////////////////////////////////////////////////////////
  // do E x \grad here also
  set_S_to_zero( S2 ) ;
  grad( S1 , S , t , mu1 ) ;
  for( i = 0 ; i < LCU ; i++ ) {
    multab( (void*)S2[i].D[0] , (void*)Fmunu[i].O[Fidx1] , (void*)S1[i].D[0] ) ;
    multab( (void*)S2[i].D[1] , (void*)Fmunu[i].O[Fidx1] , (void*)S1[i].D[1] ) ;
    multab( (void*)S2[i].D[2] , (void*)Fmunu[i].O[Fidx1] , (void*)S1[i].D[2] ) ;
    multab( (void*)S2[i].D[3] , (void*)Fmunu[i].O[Fidx1] , (void*)S1[i].D[3] ) ;
  }
  // multiply by sigma
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      S3[i].D[0][j] += mult_map[0] * S2[i].D[ sigma_map[0] ][j] ;
      S3[i].D[1][j] += mult_map[1] * S2[i].D[ sigma_map[1] ][j] ;
      S3[i].D[2][j] += mult_map[2] * S2[i].D[ sigma_map[2] ][j] ;
      S3[i].D[3][j] += mult_map[3] * S2[i].D[ sigma_map[3] ][j] ;
    }
  }
  // and the minus term
  set_S_to_zero( S2 ) ;
  grad( S1 , S , t , mu2 ) ;
  for( i = 0 ; i < LCU ; i++ ) {
    multab( (void*)S2[i].D[0] , (void*)Fmunu[i].O[Fidx2] , (void*)S1[i].D[0] ) ;
    multab( (void*)S2[i].D[1] , (void*)Fmunu[i].O[Fidx2] , (void*)S1[i].D[1] ) ;
    multab( (void*)S2[i].D[2] , (void*)Fmunu[i].O[Fidx2] , (void*)S1[i].D[2] ) ;
    multab( (void*)S2[i].D[3] , (void*)Fmunu[i].O[Fidx2] , (void*)S1[i].D[3] ) ;
  }
  // multiply by sigma
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      S3[i].D[0][j] += mult_map[0] * S2[i].D[ sigma_map[0] ][j] ;
      S3[i].D[1][j] += mult_map[1] * S2[i].D[ sigma_map[1] ][j] ;
      S3[i].D[2][j] += mult_map[2] * S2[i].D[ sigma_map[2] ][j] ;
      S3[i].D[3][j] += mult_map[3] * S2[i].D[ sigma_map[3] ][j] ;
    }
  }
  return ;
}

static void
sigmaB_G( struct halfspinor *H ,
	  struct halfspinor *S1 ,
	  const struct field Fmunu ,
	  const struct halfspinor S ,
	  const double fac )
{
  double complex A[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex sum[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
    
  // initialise sigma.B into S1
  size_t d1 , d2 , d3 , j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    S1 -> D[0][j] =  Fmunu.O[2][j] ;
    S1 -> D[1][j] =  Fmunu.O[0][j] - I * Fmunu.O[1][j] ;
    S1 -> D[2][j] =  Fmunu.O[0][j] + I * Fmunu.O[1][j] ;
    S1 -> D[3][j] = -Fmunu.O[2][j] ;
  }

  // multiply this guy with F -> S
  for( d1 = 0 ; d1 < 2 ; d1++ ) {
    for( d2 = 0 ; d2 < 2 ; d2++ ) {
	
      for( d3 = 0 ; d3 < NCNC ; d3++ ) { sum[ d3 ] = 0.0 ; }

      multab( (void*)A , (void*)S1 -> D[ d1 ] , (void*)S.D[ 2*d2 ] ) ;
      add_mat( (void*)sum , (void*)A ) ;

      multab( (void*)A , (void*)S1 -> D[ d1 + 2 ] , (void*)S.D[ 1 + 2*d2 ] ) ;
      add_mat( (void*)sum , (void*)A ) ;
      
      // update the hamiltonian
      for( j = 0 ; j < NCNC ; j++ ) {
	H -> D[ d2 + d1*2 ][j] += fac * sum[ j ] ;
      }
    }
  }
  
  return ;
}

// result is in F -> S3
static void
sigma_gradxE( struct halfspinor *S3 ,
	      struct halfspinor *S2 ,
	      struct halfspinor *S1 ,
	      struct field *Fmunu ,
	      const struct halfspinor *S ,
	      const size_t t )
{
  set_S_to_zero( S3 ) ;

  // does \sigma_x ( \grad_1 F_23 - \grad_2 F_13 )
  const size_t sigma_x[NS] = { 2 , 3 , 0 , 1 } ;
  const double complex multx[NS] = { 1 , 1 , 1 , 1 } ;
  sigma_dot_grad_x_E( S3 , S2 , S1 , Fmunu , S ,
		      t , 5 , 4 , 1 , 2 , sigma_x , multx ) ;
  
  // does \sigma_y ( \grad_2 F_03 - \grad_0 F_23 )
  const size_t sigma_y[NS] = { 2 , 3 , 0 , 1 } ;
  const double complex multy[ND] = { -I ,-I ,I ,I } ;
  sigma_dot_grad_x_E(  S3 , S2 , S1 , Fmunu , S ,
		       t , 3 , 5 , 2 , 0 , sigma_y , multy ) ;
  
  // does \sigma_z ( \grad_0 F_13 - \grad_1 F_03 )
  const size_t sigma_z[NS] = { 0 , 1 , 2 , 3} ;
  const double complex multz[ND] = { 1 , 1 ,-1 ,-1 } ;
  sigma_dot_grad_x_E( S3 , S2 , S1 , Fmunu , S ,
		      t , 4 , 3 , 0 , 1 , sigma_z , multz ) ;
  return ;
}

// this one seems tricky \sigma.( E x \grad - \grad x E ) G
// important to note that the indices have incorporated the implicit minus sign!
void
term_C3( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C3 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C3 / ( 2. * pow( 2*NRQCD.M_0 , 2 ) ) ;

  sigma_gradxE( F -> S3 , F -> S2 , F -> S1 , F -> Fmunu , F -> S , t ) ;

  update_H( F -> H , fac , F -> S3 ) ;
  
  return ;
}

// sigma.B*G
void
term_C4( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C4 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C4 / ( 2. * NRQCD.M_0 ) ;

  size_t i ;
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
  if( fabs( NRQCD.C7 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C7 / ( pow( 2*NRQCD.M_0 , 3 ) ) ;

  size_t i , mu ;
  // first term \grad^2(sigma.B.G)
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S1 , F -> S , t , mu ) ;
  }
  for ( i = 0  ; i < LCU ; i++ ) {
    sigmaB_G( &F -> H[i] , &F -> S2[i] , F -> Fmunu[i] , F -> S1[i] , fac ) ;
  }

  // second term sigma.B.\grad^2(G)
  set_S_to_zero( F -> S1 ) ;
  set_S_to_zero( F -> S2 ) ;
  for ( i = 0  ; i < LCU ; i++ ) {
    sigmaB_G( &F -> S2[i] , &F -> S1[i] , F -> Fmunu[i] , F -> S[i] , fac ) ;
  }
  set_S_to_zero( F -> S2 ) ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S2 , F -> S1 , t , mu ) ;
  }
  update_H( F -> H , fac , F -> S2 ) ;
  
  return ;
}

// this term is even worse than c3 -> use yet another temporary
void
term_C8( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C8 ) < 1E-12 ) return ;
  
  const double fac = -3.0 * NRQCD.C8 / ( 4. * pow( 2*NRQCD.M_0 , 4 ) ) ;

  size_t mu ;
  sigma_gradxE( F -> S3 , F -> S2 , F -> S1 , F -> Fmunu , F -> S , t ) ;
  
  set_S_to_zero( F -> S4 ) ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S4 , F -> S3 , t , mu ) ;
  }
  update_H( F -> H , fac , F -> S4 ) ;

  set_S_to_zero( F -> S3 ) ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S3 , F -> S , t , mu ) ;
  }
  sigma_gradxE( F -> S4 , F -> S2 , F -> S1 , F -> Fmunu , F -> S3 , t ) ;

  update_H( F -> H , fac , F -> S4 ) ;
  
  return ;
}
