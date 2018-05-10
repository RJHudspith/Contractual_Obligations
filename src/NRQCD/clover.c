/**
   @file clover.c
   @brief compute the clover field strength tensor
 */
#include "common.h"

#include "matrix_ops.h"

// the hermitian projection of the logarithm
static void
Hermitian_proj( double complex Q[ NCNC ] , 
		const double complex U[ NCNC ] )
{
#if NC == 3
  register const double cimU0 = cimag( *( U + 0 ) ) ;
  register const double cimU4 = cimag( *( U + 4 ) ) ;
  register const double cimU8 = cimag( *( U + 8 ) ) ;
  *( Q + 0 ) = ( 2 * cimU0 - cimU4 - cimU8 ) / 3. ; 
  *( Q + 1 ) = -I * ( U[1] - conj( U[3] ) ) / 2. ;  
  *( Q + 2 ) = -I * ( U[2] - conj( U[6] ) ) / 2. ; 
  *( Q + 3 ) = conj( Q[1] ) ; 
  *( Q + 4 ) = ( 2 * cimU4 - cimU0 - cimU8 ) / 3. ;
  *( Q + 5 ) = -I * ( U[5] - conj( U[7] ) ) / 2. ; 
  *( Q + 6 ) = conj( Q[2] ) ;  
  *( Q + 7 ) = conj( Q[5] ) ;  
  *( Q + 8 ) = -creal( Q[0] ) - creal( Q[4] ) ; 
#elif NC == 2
  *( Q + 0 ) = cimag( U[0] ) ;
  *( Q + 1 ) = -I * U[1] ; //OneOI2 * ( U[1] - conj( U[2] ) ) ;  
  *( Q + 2 ) = conj( Q[1] )  ; 
  *( Q + 3 ) = -Q[0] ;  
#else
  int i , j ;
  register double tr = 0.0 ;
  for( i = 0 ; i < NC ; i++ ) {
    Q[ i*(NC+1) ] = cimag( U[ i*(NC+1) ] ) ;
    tr += creal( Q[ i*(NC+1) ] ) ; 
    for( j = i+1 ; j < NC ; j++ ) {
      Q[ j+i*NC ] = -I * ( U[ j+i*NC ] - conj( U[ i+j*NC ] ) ) / 2. ;
      Q[ i+j*NC ] = conj( Q[ j+i*NC ] ) ;
    }
  }
  tr /= NC ;
  for( i = 0 ; i < NC ; i++ ) {
    Q[ i*(NC+1) ] -= tr ;
  }
#endif
  return ;
}

static void
compute_clover( double complex *F ,
		const struct site *lat ,
		const size_t i ,
		const size_t mu ,
		const size_t nu )
{
  size_t j ;
  double complex u[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex v[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  double complex sum[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  size_t s1 , s2 , s3 ;

  zero_colormatrix( sum ) ;
  
  // top right
  s1 = lat[i].neighbor[mu] ;
  s2 = lat[i].neighbor[nu] ;
  multab( (void*)u , (void*)lat[i].O[mu] , (void*)lat[s1].O[nu] ) ;
  multab_dag( (void*)v , (void*)u , (void*)lat[s2].O[mu] ) ;
  multab_dag( (void*)u , (void*)v , (void*)lat[i].O[nu] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    sum[j] += u[j] ;
  }
  // bottom right
  s1 = lat[i].back[nu] ;
  s2 = lat[s1].neighbor[mu] ;
  multabdag( (void*)u , (void*)lat[s1].O[nu] , (void*)lat[s1].O[mu] ) ;
  multab( (void*)v , (void*)u , (void*)lat[s2].O[nu] ) ;
  multab_dag( (void*)u , (void*)v , (void*)lat[i].O[mu] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    sum[j] += u[j] ;
  }
  // bottom left
  s1 = lat[i].back[mu] ;
  s2 = lat[s1].back[nu] ;
  s3 = lat[i].back[nu] ;
  multab_dagdag( (void*)u , (void*)lat[s1].O[mu] , (void*)lat[s2].O[nu] ) ;
  multab( (void*)v , (void*)u , (void*)lat[s2].O[mu] ) ;
  multab( (void*)u , (void*)v , (void*)lat[s3].O[nu] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    sum[j] += u[j] ;
  }
  // top left
  s2 = lat[i].back[mu] ;
  s1 = lat[s2].neighbor[nu] ;
  multab_dag( (void*)u , (void*)lat[i].O[nu] , (void*)lat[s1].O[mu] ) ;
  multab_dag( (void*)v , (void*)u , (void*)lat[s2].O[nu] ) ;
  multab( (void*)u , (void*)v , (void*)lat[s2].O[mu] ) ;
  for( j = 0 ; j < NCNC ; j++ ) {
    sum[j] += u[j] ;
    // clover norm, apparently a convention is to give this a minus sign
    sum[j] = -sum[j] / 4. ;
  }
  // compute traceless F - F^\dagger
  Hermitian_proj( F , sum ) ;
  return ;
}

static void
improve_cloverB( double complex *res , 
		 const struct field *F ,
		 const struct site *lat ,
		 const size_t i ,
		 const size_t t ,
		 const size_t mu ,
		 const size_t nu )
{
  double complex sum[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp1[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp2[ NCNC ] __attribute__((aligned(ALIGNMENT)));

  zero_colormatrix( sum ) ;
  
  const size_t idx = i + LCU*t ;       // is the global idx
  size_t fwd  = lat[i].neighbor[mu] ; // local Fmunu idx
  size_t bck  = lat[i].back[mu] ;
  size_t bck2 = lat[idx].back[mu] ;
  
  multab_dag( (void*)temp1 , (void*)F[fwd].O[mu] , (void*)lat[idx].O[mu] ) ;
  multab( (void*)temp2 , (void*)lat[idx].O[mu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)F[bck].O[mu] , (void*)lat[bck2].O[mu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck2].O[mu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  // there is a trick here as it wants F_{nu,mu} but this is
  // just the same as -F_{mu,nu}!
  fwd  = lat[i].neighbor[nu] ; // local Fmunu idx
  bck  = lat[i].back[nu] ;
  bck2 = lat[idx].back[nu] ;
  
  multab_dag( (void*)temp1 , (void*)F[fwd].O[mu] , (void*)lat[idx].O[nu] ) ;
  multab( (void*)temp2 , (void*)lat[idx].O[nu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)F[bck].O[mu] , (void*)lat[bck2].O[nu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck2].O[nu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    res[j] = 5*F[i].O[mu][j]/3. - sum[j]/6. ;
  }
  
  return ;
}

static void
improve_cloverE( double complex *res ,
		 struct field *F ,
		 const struct site *lat ,
		 const size_t i ,
		 const size_t t , 
		 const size_t mu ,
		 const size_t nu )
{
  double complex sum[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp1[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp2[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex Eup[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex Edn[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  
  zero_colormatrix( sum ) ;

  // mu index is always on this time plane
  const size_t idx = i + LCU*t ;       // is the global idx
  size_t fwd  = lat[i].neighbor[mu] ; // local Fmunu idx
  size_t bck  = lat[i].back[mu] ;
  size_t bck2 = lat[idx].back[mu] ;
  
  multab_dag( (void*)temp1 , (void*)F[fwd].O[mu] , (void*)lat[idx].O[mu] ) ;
  multab( (void*)temp2 , (void*)lat[idx].O[mu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)F[bck].O[mu] , (void*)lat[bck2].O[mu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck2].O[mu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  // nu index is in the one above, so we just compute these
  compute_clover( Eup , lat , lat[idx].neighbor[nu] , mu , nu ) ;
  compute_clover( Edn , lat , lat[idx].back[nu] , mu , nu ) ;

  // indexing in the nu direction
  bck  = lat[i].back[nu] ;
  bck2 = lat[idx].back[nu] ;
  
  multab_dag( (void*)temp1 , (void*)Eup , (void*)lat[idx].O[nu] ) ;
  multab( (void*)temp2 , (void*)lat[idx].O[nu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)Edn , (void*)lat[bck2].O[nu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck2].O[nu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  // improve the E-fields
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    res[j] = 5*F[i].O[mu][j]/3. - sum[j]/6. ;
  }
  
  return ;
}

// computes the E and B fields needed in the NRQCD action
void
compute_clovers( struct NRQCD_fields *F ,
		 const struct site *lat ,
		 const size_t t )
{
  const size_t idx = LCU*t ;
  size_t i ;
  
  // precompute field strength tensor if this is really expensive
  // we can do it elsewhere
  for( i = 0 ; i < LCU ; i++ ) {
    // B fields
    compute_clover( F -> Fmunu[i].O[0] , lat , i + idx , 0 , 1 ) ;
    compute_clover( F -> Fmunu[i].O[1] , lat , i + idx , 1 , 2 ) ;
    compute_clover( F -> Fmunu[i].O[2] , lat , i + idx , 2 , 0 ) ;
    // E fields
    compute_clover( F -> Fmunu[i].O[3] , lat , i + idx , 0 , 3 ) ;
    compute_clover( F -> Fmunu[i].O[4] , lat , i + idx , 1 , 3 ) ;
    compute_clover( F -> Fmunu[i].O[5] , lat , i + idx , 2 , 3 ) ;
	
  }

  // do the clover field improvement, put it into S3 and S4 as they
  // are sinks for temporaries anyway
  for( i = 0 ; i < LCU ; i++ ) {
    // bfields are easier
    improve_cloverB( F -> S3[i].D[0] , F -> Fmunu , lat , i , t , 0 , 1 ) ;
    improve_cloverB( F -> S3[i].D[1] , F -> Fmunu , lat , i , t , 1 , 2 ) ;
    improve_cloverB( F -> S3[i].D[2] , F -> Fmunu , lat , i , t , 2 , 0 ) ;
    // Efields
    improve_cloverE( F -> S4[i].D[0] , F -> Fmunu , lat , i , t , 0 , 3 ) ;
    improve_cloverE( F -> S4[i].D[1] , F -> Fmunu , lat , i , t , 1 , 3 ) ;
    improve_cloverE( F -> S4[i].D[2] , F -> Fmunu , lat , i , t , 2 , 3 ) ;
  }

  // copy back into Fmunu
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_equiv( (void*)F -> Fmunu[i].O[0] , (void*)F -> S3[i].D[0] ) ;
    colormatrix_equiv( (void*)F -> Fmunu[i].O[1] , (void*)F -> S3[i].D[1] ) ;
    colormatrix_equiv( (void*)F -> Fmunu[i].O[2] , (void*)F -> S3[i].D[2] ) ;
    // E fields
    colormatrix_equiv( (void*)F -> Fmunu[i].O[3] , (void*)F -> S4[i].D[0] ) ;
    colormatrix_equiv( (void*)F -> Fmunu[i].O[4] , (void*)F -> S4[i].D[1] ) ;
    colormatrix_equiv( (void*)F -> Fmunu[i].O[5] , (void*)F -> S4[i].D[2] ) ;
  }

  // do a test
#if 0
  double complex temp1[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp2[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex sum[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  zero_colormatrix( sum ) ;
  
  size_t mu , nu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = mu+1 ; nu < ND ; nu++ ) {
      compute_clover( temp1 , lat , i , mu , nu ) ;
      compute_clover( temp2 , lat , i , nu , mu ) ;
      add_mat( (void*)sum , temp1 ) ;
      add_mat( (void*)sum , temp2 ) ;
      print_colormatrix( temp1 ) ;
      print_colormatrix( temp2 ) ;
      print_colormatrix( sum ) ;
    }
  }
#endif

  return ;
}


