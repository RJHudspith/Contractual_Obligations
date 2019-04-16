/**
   @file clover.c
   @brief compute the clover field strength tensor
 */
#include "common.h"

#include "matrix_ops.h"  // add_mat
#include "mmul.h"        // multab

// the hermitian projection of the logarithm
static void
Hermitian_proj( double complex Q[ NCNC ] , 
		const double complex U[ NCNC ] )
{
#if NC == 3
  register const double cimU0 = cimag( *( U + 0 ) ) ;
  register const double cimU4 = cimag( *( U + 4 ) ) ;
  register const double cimU8 = cimag( *( U + 8 ) ) ;
  #ifdef HAVE_IMMINTRIN_H
  const __m128d *pU = (const __m128d*)U ;
  __m128d *pQ = (__m128d*)Q ;
  const __m128d half = _mm_set_pd( 0.5 , 0.5 ) ;
  pQ[0] = _mm_setr_pd( ( 2. * cimU0 - cimU4 - cimU8 )/3. , 0 ) ;
  pQ[1] = SSE2_iMUL( _mm_mul_pd( half ,
				 _mm_sub_pd( SSE2_CONJ( pU[3] ) ,
					     pU[1] ) ) ) ;
  pQ[2] = SSE2_iMUL( _mm_mul_pd( half ,
				  _mm_sub_pd( SSE2_CONJ( pU[6] ) ,
					      pU[2] ) ) ) ;
  pQ[3] = SSE2_iMUL( _mm_mul_pd( half ,
				 _mm_sub_pd( SSE2_CONJ( pU[1] ) ,
					     pU[3] ) ) ) ;
  pQ[4] = _mm_setr_pd( ( 2. * cimU4 - cimU0 - cimU8 )/3. , 0 ) ;
  pQ[5] = SSE2_iMUL( _mm_mul_pd( half ,
				 _mm_sub_pd( SSE2_CONJ( pU[7] ) ,
					     pU[5] ) ) ) ;
  pQ[6] = SSE2_iMUL( _mm_mul_pd( half ,
				 _mm_sub_pd( SSE2_CONJ( pU[2] ) ,
					     pU[6] ) ) ) ;
  pQ[7] = SSE2_iMUL( _mm_mul_pd( half ,
				 _mm_sub_pd( SSE2_CONJ( pU[5] ) ,
					     pU[7] ) ) ) ;
  pQ[8] = _mm_setr_pd( ( 2*cimU8 - cimU4 - cimU0 )/3. , 0 ) ;
  #else
  *( Q + 0 ) = ( 2. * cimU0 - cimU4 - cimU8 ) / 3. ;
  *( Q + 1 ) = I * 0.5 * ( conj( U[3] ) - U[1] ) ;  
  *( Q + 2 ) = I * 0.5 * ( conj( U[6] ) - U[2] ) ;
  *( Q + 3 ) = I * 0.5 * ( conj( U[1] ) - U[3] ) ;
  *( Q + 4 ) = ( 2. * cimU4 - cimU0 - cimU8 ) / 3. ;
  *( Q + 5 ) = I * 0.5 * ( conj( U[7] ) - U[5] ) ;
  *( Q + 6 ) = I * 0.5 * ( conj( U[2] ) - U[6] ) ;  
  *( Q + 7 ) = I * 0.5 * ( conj( U[5] ) - U[7] ) ;
  *( Q + 8 ) = -*( Q + 0 ) - *( Q + 4 ) ;
  #endif 
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

// does what it says and computes the clover
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
  
  // top right
  s1 = lat[i].neighbor[mu] ;
  s2 = lat[i].neighbor[nu] ;
  multab( (void*)u , (void*)lat[i].O[mu] , (void*)lat[s1].O[nu] ) ;
  multab_dag( (void*)v , (void*)u , (void*)lat[s2].O[mu] ) ;
  multab_dag( (void*)sum , (void*)v , (void*)lat[i].O[nu] ) ;

  // bottom right
  s1 = lat[i].back[nu] ;
  s2 = lat[s1].neighbor[mu] ;
  multabdag( (void*)u , (void*)lat[s1].O[nu] , (void*)lat[s1].O[mu] ) ;
  multab( (void*)v , (void*)u , (void*)lat[s2].O[nu] ) ;
  multab_dag( (void*)u , (void*)v , (void*)lat[i].O[mu] ) ;
  add_mat( (void*)sum , (void*)u ) ;

  // bottom left
  s1 = lat[i].back[mu] ;
  s2 = lat[s1].back[nu] ;
  s3 = lat[i].back[nu] ;
  multab_dagdag( (void*)u , (void*)lat[s1].O[mu] , (void*)lat[s2].O[nu] ) ;
  multab( (void*)v , (void*)u , (void*)lat[s2].O[mu] ) ;
  multab( (void*)u , (void*)v , (void*)lat[s3].O[nu] ) ;
  add_mat( (void*)sum , (void*)u ) ;

  // top left
  s2 = lat[i].back[mu] ;
  s1 = lat[s2].neighbor[nu] ;
  multab_dag( (void*)u , (void*)lat[i].O[nu] , (void*)lat[s1].O[mu] ) ;
  multab_dag( (void*)v , (void*)u , (void*)lat[s2].O[nu] ) ;
  multab( (void*)u , (void*)v , (void*)lat[s2].O[mu] ) ;
  add_mat( (void*)sum , (void*)u ) ;

  // clover norm, apparently a convention is to give this a minus sign
  for( j = 0 ; j < NCNC ; j++ ) {
    sum[j] = -sum[j] / 4. ;
  }
  
  // compute traceless F - F^\dagger
  Hermitian_proj( F , sum ) ;
  
  return ;
}

// computes the improved clover of lepage and magnea
static void
improve_clover( double complex *res ,
		const struct site *lat ,
		const size_t i ,
		const size_t mu ,
		const size_t nu ,
		const double efac )
{  
  double complex sum[ NCNC ]   __attribute__((aligned(ALIGNMENT)));
  double complex temp1[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex temp2[ NCNC ] __attribute__((aligned(ALIGNMENT)));
  double complex Eup[ NCNC ]   __attribute__((aligned(ALIGNMENT)));
  double complex Edn[ NCNC ]   __attribute__((aligned(ALIGNMENT)));
  
  // standard clover term
  compute_clover( res , lat , i , mu , nu ) ;
  
  zero_colormatrix( sum ) ;

  // mu index contributions
  const size_t bck_mu = lat[i].back[mu] ;
  compute_clover( Eup , lat , lat[i].neighbor[mu] , mu , nu ) ;
  compute_clover( Edn , lat , lat[i].back[mu] , mu , nu ) ;
  
  multab_dag( (void*)temp1 , (void*)Eup , (void*)lat[i].O[mu] ) ;
  multab( (void*)temp2 , (void*)lat[i].O[mu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)Edn , (void*)lat[bck_mu].O[mu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck_mu].O[mu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  // now do the nu index
  const size_t bck_nu = lat[i].back[nu] ;
  compute_clover( Eup , lat , lat[i].neighbor[nu] , mu , nu ) ;
  compute_clover( Edn , lat , bck_nu , mu , nu ) ;
  
  multab_dag( (void*)temp1 , (void*)Eup , (void*)lat[i].O[nu] ) ;
  multab( (void*)temp2 , (void*)lat[i].O[nu] , (void*)temp1 ) ;
  add_mat( (void*)sum , (void*)temp2 ) ;
  
  multab( (void*)temp1 , (void*)Edn , (void*)lat[bck_nu].O[nu] ) ;
  multabdag( (void*)temp2 , (void*)lat[bck_nu].O[nu] , (void*)temp1 ) ; 
  add_mat( (void*)sum , (void*)temp2 ) ;

  // improve the clover fields
  size_t j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    res[j] = efac*res[j] - sum[j]/6. ;
  }
  return ;
}

// computes the E and B fields needed in the NRQCD action
void
compute_clovers( struct NRQCD_fields *F ,
		 const struct site *lat ,
		 const size_t t ,
		 const double U_0 )
{
  // e-clover factor, second term is to correct for the non-unitary ness
  const double efac = ( 4. + 1./( U_0 * U_0 ) ) / 3. ;
  
  const size_t idx = LCU*t ;
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    // B fields are defined as B_{i} = \epsilon_{ijk} F_{jk}
    improve_clover( F -> Fmunu[i].O[0] , lat , i + idx , 1 , 2 , efac ) ;
    improve_clover( F -> Fmunu[i].O[1] , lat , i + idx , 2 , 0 , efac ) ;
    improve_clover( F -> Fmunu[i].O[2] , lat , i + idx , 0 , 1 , efac ) ;
    // E fields are defined as E_{i} = F_{t i}
    improve_clover( F -> Fmunu[i].O[3] , lat , i + idx , 3 , 0 , efac ) ;
    improve_clover( F -> Fmunu[i].O[4] , lat , i + idx , 3 , 1 , efac ) ;
    improve_clover( F -> Fmunu[i].O[5] , lat , i + idx , 3 , 2 , efac ) ;
    // these last ones are used in the improved derivative
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t Uidx = i + t*LCU ;
      const size_t Ufwd = lat[ Uidx ].neighbor[mu] ;
      const size_t Ubck = lat[ Uidx ].back[mu] ;
      const size_t Ubck2 = lat[ Ubck ].back[mu] ;
      multab( (void*)F -> Fmunu[i].O[6+2*mu] , (void*)lat[ Uidx ].O[mu] , (void*)lat[ Ufwd ].O[mu] ) ;
      multab( (void*)F -> Fmunu[i].O[7+2*mu] , (void*)lat[ Ubck2 ].O[mu] , (void*)lat[ Ubck ].O[mu] ) ;
    }
  }
  
  return ;
}


