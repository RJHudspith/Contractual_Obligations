/**
   @file spin_independent.c
   @brief collection of terms in the NRQCD hamiltonian
 */
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

// This part computes -C_0 * ( grad^2 / 2*M_0 ) S(x,t) 
void
term_C0( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD )
{
  if( fabs( NRQCD.C0 ) < 1E-12 ) return ;
  
  const double fac = -NRQCD.C0 / ( 2. * NRQCD.M_0 ) ;

  set_S_to_zero( F -> S1 ) ;

  grad_sq( F -> S1 , F -> S , t ) ;
  
  update_H( F -> H , fac , F -> S1 ) ;

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

  const double unicf = 0.0 ;
  
  set_S_to_zero( F -> S1 ) ;
  set_S_to_zero( F -> S2 ) ;

  size_t i , mu , j ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S1 , F -> S , t , mu ) ;
  }

  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S2  , F -> S1 , t , mu ) ;
  }

  // multiply by factor
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      F -> H[i].D[0][j] += fac * ( F -> S2[i].D[0][j] + unicf * F->S[i].D[0][j] ) ;
      F -> H[i].D[1][j] += fac * ( F -> S2[i].D[1][j] + unicf * F->S[i].D[1][j] ) ;
      F -> H[i].D[2][j] += fac * ( F -> S2[i].D[2][j] + unicf * F->S[i].D[2][j] ) ;
      F -> H[i].D[3][j] += fac * ( F -> S2[i].D[3][j] + unicf * F->S[i].D[3][j] ) ;
    }
  }
    
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

  set_S_to_zero( F -> S1 ) ;
  set_S_to_zero( F -> S2 ) ;
  
  size_t mu , i , j ;
  // compute i.Delta.E - E.Delta
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    // E * prop
    for( i = 0 ; i < LCU ; i++ ) {
      multab( (void*)F -> S1[i].D[0] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] , (void*)F -> S[i].D[0] ) ;
      multab( (void*)F -> S1[i].D[1] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] , (void*)F -> S[i].D[1] ) ;
      multab( (void*)F -> S1[i].D[2] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] , (void*)F -> S[i].D[2] ) ;
      multab( (void*)F -> S1[i].D[3] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] , (void*)F -> S[i].D[3] ) ;
    }
    // call derivative into S2
    grad( F -> S2 , F -> S1 , t , mu ) ;
    
    // update the Hamiltonian
    for( i = 0 ; i < LCU ; i++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	F -> H[i].D[0][j] += I * fac * F -> S2[i].D[0][j] ;
	F -> H[i].D[1][j] += I * fac * F -> S2[i].D[1][j] ;
	F -> H[i].D[2][j] += I * fac * F -> S2[i].D[2][j] ;
	F -> H[i].D[3][j] += I * fac * F -> S2[i].D[3][j] ;
      }
    }
  }

  set_S_to_zero( F -> S1 ) ;
  
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    // call improved derivative into S1
    grad( F -> S1 , F -> S , t , mu ) ;
    // E * prop
    for( i = 0 ; i < LCU ; i++ ) {
      multab( (void*)F -> S2[i].D[0] , (void*)F -> S[i].D[0] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] ) ;
      multab( (void*)F -> S2[i].D[1] , (void*)F -> S[i].D[1] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] ) ;
      multab( (void*)F -> S2[i].D[2] , (void*)F -> S[i].D[2] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] ) ;
      multab( (void*)F -> S2[i].D[3] , (void*)F -> S[i].D[3] , (void*)F -> Fmunu[i].O[ (ND-1) + mu ] ) ;  
    }
    // update the Hamiltonian
    for( i = 0 ; i < LCU ; i++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	F -> H[i].D[0][j] -= I * fac * F -> S2[i].D[0][j] ;
	F -> H[i].D[1][j] -= I * fac * F -> S2[i].D[1][j] ;
	F -> H[i].D[2][j] -= I * fac * F -> S2[i].D[2][j] ;
	F -> H[i].D[3][j] -= I * fac * F -> S2[i].D[3][j] ;
      }
    }
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

  const double unicf = 0.0 ;

  set_S_to_zero( F -> S1 ) ;
  set_S_to_zero( F -> S2 ) ;
  
  size_t i , mu , j ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    grad2( F -> S1 , F -> S , t , mu ) ;
    grad2( F -> S2 , F -> S1 , t , mu ) ;
  }

  // multiply by factor
  for( i = 0 ; i < LCU ; i++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      F -> H[i].D[0][j] += fac * ( F -> S2[i].D[0][j] + unicf * F->S[i].D[0][j] ) ;
      F -> H[i].D[1][j] += fac * ( F -> S2[i].D[1][j] + unicf * F->S[i].D[1][j] ) ;
      F -> H[i].D[2][j] += fac * ( F -> S2[i].D[2][j] + unicf * F->S[i].D[2][j] ) ;
      F -> H[i].D[3][j] += fac * ( F -> S2[i].D[3][j] + unicf * F->S[i].D[3][j] ) ;
    }
  }
  
  return ;
}
