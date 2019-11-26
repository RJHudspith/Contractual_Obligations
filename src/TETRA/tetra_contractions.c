/**
   @file tetra_contractions.c
   @brief tetraquark contractions
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_r()
#include "gammas.h"             // Cgmu()
#include "spinmatrix_ops.h"     // trace_prod, get_spinmatrix()
#include "spinor_ops.h"         // transpose_spinor()
#include "tetra_contractions.h" // alphabetising
#include "penta_contractions.h" // idx()
#include "Ospinor.h"            // spinor_to_Ospinor()

// block of possible gammas in the contraction
struct tblock { 
  struct gamma G5 ;    // light diquark gamma
  struct gamma Gi ;    // heavy diquark gamma
  struct gamma t_G5 ;  // gt ( G5^\dagger ) gt
  struct gamma t_Gi ;  // gt ( Gi^\dagger ) gt
  struct gamma CG5 ;   // C.G5
  struct gamma CGi ;   // C.Gi
  struct gamma t_CG5 ; // gt ( (C.G5)^\dagger ) gt
  struct gamma t_CGi ; // gt ( (C.Gi)^\dagger ) gt
} ;

// perform the contraction
static double complex
contract_O1O1( const struct block *C1 , 
	       const struct block *C2 ,
	       const GLU_bool H1H2_degenerate )
{
  register double complex sum = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;

    // Hsum is the sum of heavy contractions
    register double complex Hsum =
      spinmatrix_trace( C2[ element( c , a , d , b ) ].M ) ;
    
    // cross color term by interchanging Bs picks up a minus
    if( H1H2_degenerate == GLU_TRUE ) {
      Hsum -= spinmatrix_trace( C2[ element( c , b , d , a ) ].M ) ;
    // second terms come from epsilon identity
    } else {
      Hsum -= spinmatrix_trace( C2[ element( d , a , c , b ) ].M ) ;
      Hsum -= spinmatrix_trace( C2[ element( c , b , d , a ) ].M ) ;
      Hsum += spinmatrix_trace( C2[ element( d , b , c , a ) ].M ) ;
    }

    // multiply by the light contraction
    sum += spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) * Hsum ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*sum : sum ;
}

// diquark-dimeson cross term
static double complex
contract_O1O2_1( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;

    // single dirac trace
    sum += 
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , a , c , b ) ].M ) ;
    // cross term
    sum -=
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , b , c , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*sum : sum ;
}

// diquark-dimeson cross term
static double complex
contract_O1O2_2( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // single dirac trace
    sum += 
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , a , d , b ) ].M ) ;
    sum -=
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , b , d , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*sum : sum ;
}

// dimeson-diquark cross term
static double complex
contract_O2O1_1( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum += 
      trace_prod_spinmatrices( C1[ element( c , b , b , d ) ].M ,
			       C2[ element( a , c , d , a ) ].M ) ;
    sum -= 
      trace_prod_spinmatrices( C1[ element( d , b , b , d ) ].M ,
			       C2[ element( a , c , c , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*sum : sum ;
}

// dimeson-diquark cross term
static double complex
contract_O2O1_2( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum += 
      trace_prod_spinmatrices( C1[ element( c , a , a , d ) ].M ,
			       C2[ element( b , c , d , b ) ].M ) ;
    // cross term
    sum -= 
      trace_prod_spinmatrices( C1[ element( d , a , a , d ) ].M ,
			       C2[ element( b , c , c , b ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*sum : sum ;
}

static double complex
contract_O2O2_1( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum1 += 
      spinmatrix_trace( C1[ element( c , a , a , c ) ].M ) *
      spinmatrix_trace( C2[ element( d , b , b , d ) ].M ) ;
    // cross term
    sum2 += 
      trace_prod_spinmatrices( C1[ element( d , a , a , c ) ].M ,
			       C2[ element( c , b , b , d ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? ( sum1 - sum2 ) : sum1 ;
}

static double complex
contract_O2O2_2( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum1 += 
      spinmatrix_trace( C1[ element( d , a , a , d ) ].M ) *
      spinmatrix_trace( C2[ element( c , b , b , c ) ].M ) ;
    // cross term
    sum2 += 
      trace_prod_spinmatrices( C1[ element( c , a , a , d ) ].M ,
			       C2[ element( d , b , b , c ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? ( sum1 - sum2 ) : sum2 ;
}

// get an element from a b c d
size_t
element( const size_t a , const size_t b , const size_t c , const size_t d ) 
{
  return d + NC * ( c + NC * ( b + NC * a ) ) ;
}

// little function to determine a,b,c and d. d runs fastest, a slowest
void
get_abcd( size_t *a , size_t *b , size_t *c , size_t *d , const size_t abcd )
{
  // cache this
  const size_t NCNCNC = NCNC * NC ;
  *a = ( abcd ) / ( NCNCNC ) ;
  *b = ( abcd % NCNCNC ) / ( NCNC ) ;
  *c = ( abcd % NCNC ) / ( NC ) ;
  *d = ( abcd % NC ) ;
}

// precompute a "block" as an array of spinmatrices
// with open colors a,b,c,d in a flattened array
void
precompute_block( struct block *C1 ,
		  const struct spinor S1 ,
		  const struct gamma G1 ,
		  const struct spinor S2 ,
		  const struct gamma G2 )
{
  struct spinor S1g1 = S1 , S2g2 = S2 ;
  gamma_mul_r( &S1g1 , G1 ) ;
  gamma_mul_r( &S2g2 , G2 ) ;
  double complex D1[ NSNS ] , D2[ NSNS ] ;
  size_t ab , cd , a , b , c , d ;
  for( ab = 0 ; ab < NCNC ; ab++ ) {
    a = ab / NC ;
    b = ab % NC ;
    get_spinmatrix( D1 , S1g1 , a , b ) ;
    for( cd = 0 ; cd < NCNC ; cd++ ) {
      c = cd / NC ;
      d = cd % NC ;
      get_spinmatrix( D2 , S2g2 , c , d ) ;
      spinmatrix_multiply( C1 -> M , D1 , D2 ) ;
      C1++ ;
    }
  }
  return ;
}

static inline void
setM( double complex *C1 ,
      const struct spinmatrix temp )
{
#ifdef __AVX__
  double *pC1 = (double*)C1 ;
  const __m256d *pt = (const __m256d*)temp.D ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
  _mm256_storeu_pd( pC1 , *pt ) ; pC1+=4 ; pt++ ;
#else
  C1[0] = temp.D[0][0] ; C1[1] = temp.D[0][1] ; C1[2] = temp.D[0][2] ; C1[3] = temp.D[0][3] ;
  C1[4] = temp.D[1][0] ; C1[5] = temp.D[1][1] ; C1[6] = temp.D[1][2] ; C1[7] = temp.D[1][3] ;
  C1[8] = temp.D[2][0] ; C1[9] = temp.D[2][1] ; C1[10] = temp.D[2][2] ; C1[11] = temp.D[2][3] ;
  C1[12] = temp.D[3][0] ; C1[13] = temp.D[3][1] ; C1[14] = temp.D[3][2] ; C1[15] = temp.D[3][3] ;
#endif
}

// precompute a "block" as an array of spinmatrices
// with open colors a,b,c,d in a flattened array
static void
precompute_block2( struct block *C1 ,
		   const struct Ospinor S1 ,
		   const struct gamma G1 ,
		   const struct Ospinor S2 ,
		   const struct gamma G2 )
{
  struct Ospinor S1g1 = S1 , S2g2 = S2 ;
  gamma_mul_r_Ospinor( &S1g1 , G1 ) ;
  gamma_mul_r_OspinorT( &S2g2 , G2 ) ;
  size_t ab , a , b ;
  struct spinmatrix temp __attribute__ ((aligned(SPINT_ALIGNMENT))) ;
  for( ab = 0 ; ab < NCNC ; ab++ ) {
    a = ab / NC ; b = ab % NC ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[0][0].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[0][1].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[0][2].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[1][0].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[1][1].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[1][2].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[2][0].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[2][1].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;
    spinmatrix_multiply_T_avx( temp.D , S1g1.C[a][b].D , S2g2.C[2][2].D ) ;
    setM( C1 -> M , temp ) ; C1++ ;    
  }
  return ;
}

// perform the contraction of the tetraquark with appropriate mixing
// has a block-symmetric structure of TETRA_NBLOCKxTETRA_NBLOCK matrices
// Top Left is the Diquark - Diquark correlator 
// Top right is the Diquark - Dimeson mixing
// Bottom left is the Dimeson - Diquark
// Bottom right is the Dimeson - Dimeson
// when we have ND==3 we turn off the negative term in the meson-meson
int
tetras( double complex *result ,
	const struct spinor L1 , 
	const struct spinor L2 ,
	const struct spinor bwdH1 ,
	const struct spinor bwdH2 ,
	const struct gamma *GAMMAS ,
	const size_t mu , 
	const GLU_bool L1L2_degenerate ,
	const GLU_bool H1H2_degenerate )
{
#if (TETRA_NBLOCK > 8) || (TETRA_NBLOCK < 1)
  printf( stderr , "[TETRA] specified TETRA_NBLOCK not usable %d" , TETRA_NBLOCK ) ;
  return FAILURE ;
#endif
  
  // timelike gamma matrix
  const struct gamma gt = GAMMAS[ GAMMA_T ] ;

  // maps for the various second-index structure
  // the first three are the 1^+, and the last index is 0^+
  const size_t numap_Vi[ 4 ]  = { GAMMA_X , GAMMA_Y , GAMMA_Z , GAMMA_5 } ;
  const size_t numap_Tit[ 4 ] = { TXT , TYT , TZT , AT } ;
  const size_t numap_Tij[ 4 ] = { TYZ , TZX , TXY , GAMMA_T } ;
  const size_t numap_Ai[ 4 ]  = { AX , AY , AZ , IDENTITY } ;

  // b-map is a dirty hack!!
#ifdef TETRA_NRQCD_HACK
  const size_t dm[ TETRA_NBLOCK ] = { 0 , 2 } ;
  const size_t mm[ TETRA_NBLOCK ] = { 0 , 1 } ;
#else
  const size_t dm[ 8 ] = { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 } ;
  const size_t mm[ 8 ] = { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 } ;
#endif
  
  // change everything to a block - structure
  // first index is the ud - part and the second is the bb-part
  struct tblock blck[ 8 ] =
    {
      { .G5 = GAMMAS[ GAMMA_5 ]  , .Gi = GAMMAS[ numap_Vi[ mu ] ] }  , // P V_mu
      { .G5 = GAMMAS[ IDENTITY ] , .Gi = GAMMAS[ numap_Ai[ mu ] ] }  , // I A_mu
      { .G5 = GAMMAS[ AT ]       , .Gi = GAMMAS[ numap_Tit[ mu ] ] } , // At Sigma_mut
      { .G5 = GAMMAS[ GAMMA_T ]  , .Gi = GAMMAS[ numap_Tij[ mu ] ] } , // V_t T_ij
      // second block have me swapping the gammas for the heavy diquark
      { .G5 = GAMMAS[ GAMMA_5 ] ,  .Gi = GAMMAS[ numap_Tit[ mu ] ] } , // P Sigma_mut
      { .G5 = GAMMAS[ IDENTITY ] , .Gi = GAMMAS[ numap_Tij[ mu ] ] } , // I T_ij
      { .G5 = GAMMAS[ AT ] ,       .Gi = GAMMAS[ numap_Vi[ mu ] ] }  , // At V_mu
      { .G5 = GAMMAS[ GAMMA_T ] ,  .Gi = GAMMAS[ numap_Ai[ mu ] ] }  , // V_t A_mu
    } ;
  
  // transposed Ospinor temporaries
  struct Ospinor OL1  = spinor_to_Ospinor( L1 ) ;
  struct Ospinor OL2  = spinor_to_Ospinor( L2 ) ;
  struct Ospinor OL1T = spinor_to_Ospinor( transpose_spinor( L1 ) ) ;

  struct Ospinor ObwdH1 = spinor_to_Ospinor( bwdH1 ) ;
  struct Ospinor ObwdH2 = spinor_to_Ospinor( bwdH2 ) ;
  struct Ospinor ObwdH2T = spinor_to_Ospinor( transpose_spinor( bwdH2 ) ) ;

  // temporaries are 
  const size_t Nco = NCNC*NCNC ;
  struct block *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , ALIGNMENT , Nco*sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C2 , ALIGNMENT , Nco*sizeof( struct block ) ) != 0 ) {
    free( C1 ) ; free( C2 ) ;
    fprintf( stderr , "[MALLOC] failed to allocate C1/C2 temporaries" ) ;
    return sqrt(-1) ;
  }

  // compute all the usual gamma structures needed				
  size_t B1 , B2 ;
  for( B1 = 0 ; B1 < 8 ; B1++ ) {
    // tilde operations
    blck[B1].t_G5 = gt_Gdag_gt( blck[B1].G5 , gt ) ;
    blck[B1].t_Gi = gt_Gdag_gt( blck[B1].Gi , gt ) ;
    // C operations
    blck[B1].CG5 = CGmu( blck[B1].G5 , GAMMAS ) ;
    blck[B1].CGi = CGmu( blck[B1].Gi , GAMMAS ) ;
    // tilde-C operations
    blck[B1].t_CG5 = gt_Gdag_gt( blck[B1].CG5 , gt ) ;
    blck[B1].t_CGi = gt_Gdag_gt( blck[B1].CGi , gt ) ;
  }
  
  // loop blocks
  for( B1 = 0 ; B1 < TETRA_NBLOCK ; B1++ ) {
    for( B2 = 0 ; B2 < TETRA_NBLOCK ; B2++ ) {

      const size_t idx1 = B2 + 2 * (size_t)TETRA_NBLOCK*B1 ;

      // Diquark-Diquark
      precompute_block2( C1 , OL1T , blck[ dm[B1] ].CG5 , OL2 , blck[ dm[B2] ].t_CG5 ) ;
      precompute_block2( C2 , ObwdH1 , blck[ dm[B1] ].CGi , ObwdH2T , blck[ dm[ B2 ]].t_CGi ) ;
      
      result[idx1] = contract_O1O1( C1 , C2 , H1H2_degenerate ) ;

      /////////////////// Diquark-AntiDiquark - Dimeson mixing terms
      const size_t idx2 = idx1 + TETRA_NBLOCK ;

      precompute_block2( C1 , OL1T , blck[dm[B1]].CG5 , OL2 , gamma_transpose( blck[mm[B2]].t_Gi ) ) ;
      precompute_block2( C2 , ObwdH1 , blck[ dm[B1]].CGi , ObwdH2T , blck[mm[B2]].t_G5 ) ;
      result[idx2]  = contract_O1O2_1( C1 , C2 , H1H2_degenerate ) ;

      if( mu < (size_t)(ND-1) ) {
	precompute_block2( C1 , OL1T , blck[ dm[B1] ].CG5 , OL2 , gamma_transpose( blck[mm[B2]].t_G5 ) ) ;
	precompute_block2( C2 , ObwdH1 , blck[dm[B1]].CGi , ObwdH2T , blck[mm[B2]].t_Gi ) ;
	result[idx2] -= contract_O1O2_2( C1 , C2 , H1H2_degenerate ) ;
      }
      
      /////////////////// Dimeson -> Diquark Anti-Diquark mixing terms
      const size_t idx3 = 2*TETRA_NBLOCK*TETRA_NBLOCK + B2 + 2 * TETRA_NBLOCK * B1 ;
      
      // O_2 O_1 -- term 1
      precompute_block2( C1 , ObwdH1 , gamma_transpose( blck[mm[B1]].Gi ) , OL2 , blck[dm[B2]].t_CG5 ) ;
      precompute_block2( C2 , OL1T , blck[mm[B1]].G5 , ObwdH2T , blck[dm[B2]].t_CGi ) ;
      result[idx3]  = contract_O2O1_1( C1 , C2 , H1H2_degenerate ) ;

      // O_2 O_1 -- term 2 has the minus sign
      if( mu < (size_t)(ND-1) ) {
	precompute_block2( C1 , ObwdH1 , gamma_transpose( blck[mm[B1]].G5 ) , OL2 , blck[dm[B2]].t_CG5 ) ;
	precompute_block2( C2 , OL1T , blck[mm[B1]].Gi , ObwdH2T , blck[dm[B2]].t_CGi ) ;
	result[idx3] -= contract_O2O1_2( C1 , C2 , H1H2_degenerate ) ;
      }
      
      ////////////////// Dimeson -> Dimeson mixing terms
      const size_t idx4 = 2*TETRA_NBLOCK*TETRA_NBLOCK + TETRA_NBLOCK + B2 + 2 * TETRA_NBLOCK * B1 ;
      
      // O_2 O_2 -- term 1 is positive 
      precompute_block2( C1 , ObwdH1 , blck[mm[B1]].G5 , OL1 , blck[mm[B2]].t_G5 ) ;
      precompute_block2( C2 , ObwdH2 , blck[mm[B1]].Gi , OL2 , blck[mm[B2]].t_Gi ) ;
      result[idx4]  = contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;

      if( mu < (size_t)(ND-1) ) {
	// O_2 O_2 -- term 2 is -( a b^\dagger )
	precompute_block2( C1 , ObwdH1 , blck[mm[B1]].G5 , OL1 , blck[mm[B2]].t_Gi ) ;
	precompute_block2( C2 , ObwdH2 , blck[mm[B1]].Gi , OL2 , blck[mm[B2]].t_G5 ) ;  
	result[idx4] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;
	
	// need to do the others where L1 and L2 are swapped, this is only 
	// a concern for the dimeson - dimeson
	if( L1L2_degenerate == GLU_FALSE ) {
	  // O2O2 -- term 3 is -( b a^\dagger )
	  precompute_block2( C1 , ObwdH1 , blck[mm[B1]].G5 , OL2 , blck[mm[B2]].t_Gi ) ;
	  precompute_block2( C2 , ObwdH2 , blck[mm[B1]].Gi , OL1 , blck[mm[B2]].t_G5 ) ;  
	  result[idx4] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;
	  
	  // O2O2 -- term 4 is ( b b^\dagger )
	  precompute_block2( C1 , ObwdH1 , blck[mm[B1]].G5 , OL2 , blck[mm[B2]].t_G5 ) ;
	  precompute_block2( C2 , ObwdH2 , blck[mm[B1]].Gi , OL1 , blck[mm[B2]].t_Gi ) ;
	  result[idx4] += contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;
	} else {
	  result[idx4] *= 2 ;
	}
	//
      }
      // mu < ND-1
    }
  }

  free( C1 ) ; free( C2 ) ;

  return SUCCESS ;
}
