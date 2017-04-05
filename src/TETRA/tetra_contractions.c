/**
   @file tetra_contractions.c
   @brief tetraquark contractions

   TODO :: non-degenerate heavy?
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_r()
#include "gammas.h"             // Cgmu()
#include "spinmatrix_ops.h"     // trace_prod, get_spinmatrix()
#include "spinor_ops.h"         // transpose_spinor()
#include "tetra_contractions.h" // alphabetising

// positive parity operators
#define POSPOS

// perform the contraction
static double complex
contract_O1O1( const struct block *C1 , 
	       const struct block *C2 ,
	       const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual forward streaming term
    sum1 += 
      spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) *
      spinmatrix_trace( C2[ element( c , a , d , b ) ].M ) ;
    // cross color term by interchanging Bs
    sum2 += 
      spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) * 
      spinmatrix_trace( C2[ element( d , a , c , b ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*( sum1 - sum2 ) : 4*sum1 ;
}

// diquark-dimeson cross term
static double complex
contract_O1O2_1( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // single dirac trace
    sum1 += 
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , a , c , b ) ].M ) ;
    // cross term
    sum2 +=
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , b , c , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*( sum1 - sum2 ) : 2*sum1 ;
}

// diquark-dimeson cross term
static double complex
contract_O1O2_2( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // single dirac trace
    sum1 += 
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , a , d , b ) ].M ) ;
    // cross term
    sum2 +=
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , b , d , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*( sum1 - sum2 ) : 2*sum1 ;
}

// dimeson-diquark cross term
static double complex
contract_O2O1_1( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum1 += 
      trace_prod_spinmatrices( C1[ element( c , b , b , d ) ].M ,
			       C2[ element( a , c , d , a ) ].M ) ;
    // cross term
    sum2 += 
      trace_prod_spinmatrices( C1[ element( d , b , b , d ) ].M ,
			       C2[ element( a , c , c , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*( sum1 - sum2 ) : 2*sum1 ;
}

// dimeson-diquark cross term
static double complex
contract_O2O1_2( const struct block *C1 , 
		 const struct block *C2 ,
		 const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // usual meson product
    sum1 += 
      trace_prod_spinmatrices( C1[ element( c , a , a , d ) ].M ,
			       C2[ element( b , c , d , b ) ].M ) ;
    // cross term
    sum2 += 
      trace_prod_spinmatrices( C1[ element( d , a , a , d ) ].M ,
			       C2[ element( b , c , c , b ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 2*( sum1 - sum2 ) : 2*sum1 ;
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
  return ( H1H2_degenerate == GLU_TRUE ) ? ( sum1 - sum2 ) : sum1 ;
}

// mixing of heavy mesons
static double complex
contract_O3( const struct block *C1 , 
	     const struct block *C2 ,
	     const GLU_bool H1H2_degenerate )
{
  register double complex sum1 = 0.0 ;
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // meson product
    sum1 += 
      trace_prod_spinmatrices( C1[ element( d , a , a , c ) ].M ,
			       C2[ element( c , b , b , d ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return sum1 ;
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

// perform the contraction
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
  // timelike gamma matrix
  const struct gamma gt = GAMMAS[ GAMMA_T ] ;

  // New!! Map of the other gammas
#ifdef POSPOS // positive-positive parity operators
  const size_t pmap[ 2 ] = { 5 , 9 } ;                  // gamma_5, A_t
  const size_t nu1_map[ ND - 1 ] = { 0 , 1 , 2 } ;      // gamma_i
  const size_t nu2_map[ ND - 1 ] = { 12 , 14 , 15 } ;   // T_it
#else // negative parity operators
  const size_t pmap[ 2 ] = { 3 , 4 } ;                  // gamma_t, I
  const size_t nu1_map[ ND - 1 ] = { 6 , 7 , 8 } ;      // A_i
  const size_t nu2_map[ ND - 1 ] = { 10 , 11 , 13 } ;   // T_ij
#endif
  
  // rename gamma matrices to match the contractions above  
  const struct gamma gi[2] = { GAMMAS[ nu1_map[mu] ] ,
			       GAMMAS[ nu2_map[mu] ] } ;
  const struct gamma tildegi[2] = { gt_Gdag_gt( gi[0] , gt ) ,
				    gt_Gdag_gt( gi[1] , gt ) } ;
    
  const struct gamma g5[2] = { GAMMAS[ pmap[0] ] ,
			       GAMMAS[ pmap[1] ] } ;
  const struct gamma tildeg5[2] = { gt_Gdag_gt( g5[0] , gt ) ,
				    gt_Gdag_gt( g5[1] , gt ) } ;

  // charge conjugation matrices
  const struct gamma Cgi[2] = { CGmu( gi[0] , GAMMAS ) ,
				CGmu( gi[1] , GAMMAS ) } ;
  const struct gamma tildeCgi[2] = { gt_Gdag_gt( Cgi[0] , gt ) ,
				     gt_Gdag_gt( Cgi[1] , gt ) } ;

  const struct gamma Cg5[2] = { CGmu( g5[0] , GAMMAS ) ,
				CGmu( g5[1] , GAMMAS ) } ;
  const struct gamma tildeCg5[2] = { gt_Gdag_gt( Cg5[0] , gt ) ,
				     gt_Gdag_gt( Cg5[1] , gt ) } ;

  // transposed spinor temporaries
  struct spinor L1T = transpose_spinor( L1 ) ;
  struct spinor bwdH2T = transpose_spinor( bwdH2 ) ;

  // temporaries are 
  const size_t Nco = NCNC*NCNC ;
  struct block *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco*sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco*sizeof( struct block ) ) != 0 ) {
    free( C1 ) ; free( C2 ) ;
    fprintf( stderr , "[MALLOC] failed to allocate C1/C2 temporaries" ) ;
    return sqrt(-1) ;
  }

  // loop 2x2 blocks
  size_t block , nu ;
  for( block = 0 ; block < 4 ; block++ ) {
    const size_t mu1 = ( block * 4 ) / 8 ;          // is the first pseudoscalar index
    const size_t mu2 = ( ( block * 4 ) / 4 ) % 2 ;  // is the second pseudoscalar index
    for( nu = 0 ; nu < 4 ; nu++ ) {

      // second indices
      const size_t nu1 = nu/2 ;
      const size_t nu2 = nu%2 ;

      // gives the result matrix local idx
      const size_t loc_idx = mu1 * 16 + mu2 * 2 + nu1 * 8 + nu2 ;

      /////////////////// Diquark-AntiDiquark - Diquark-AntiDiquark mixing terms
      const size_t idx1 = loc_idx ;
      
      precompute_block( C1 , L1T , Cg5[mu1] , L2 , tildeCg5[mu2] ) ;
      precompute_block( C2 , bwdH1 , Cgi[nu1] , bwdH2T , tildeCgi[nu2] ) ;
      
      result[idx1] = contract_O1O1( C1 , C2 , H1H2_degenerate ) ;

      /////////////////// Diquark-AntiDiquark - Dimeson mixing terms
      const size_t idx2 = 4 + loc_idx ;
      
      precompute_block( C1 , L1T , Cg5[mu1] , L2 , gamma_transpose( tildegi[nu2] ) ) ;
      precompute_block( C2 , bwdH1 , Cgi[nu1] , bwdH2T , tildeg5[mu2] ) ;
      result[idx2]  = contract_O1O2_1( C1 , C2 , H1H2_degenerate ) ;

      precompute_block( C1 , L1T , Cg5[mu1] , L2 , gamma_transpose( tildeg5[mu2] ) ) ;
      precompute_block( C2 , bwdH1 , Cgi[nu1] , bwdH2T , tildegi[nu2] ) ;
      result[idx2] -= contract_O1O2_2( C1 , C2 , H1H2_degenerate ) ;

      /////////////////// Dimeson -> Diquark Anti-Diquark mixing terms
      const size_t idx3 = 32 + loc_idx ;
      
      // O_2 O_1 -- term 1
      precompute_block( C1 , bwdH1 , gamma_transpose( gi[nu1] ) , L2 , tildeCg5[mu2] ) ;
      precompute_block( C2 , L1T , g5[mu1] , bwdH2T , tildeCgi[nu2] ) ;
      result[idx3]  = contract_O2O1_1( C1 , C2 , H1H2_degenerate ) ;

      // O_2 O_1 -- term 2 has the minus sign
      precompute_block( C1 , bwdH1 , gamma_transpose( g5[mu1] ) , L2 , tildeCg5[mu2] ) ;
      precompute_block( C2 , L1T , gi[nu1] , bwdH2T , tildeCgi[nu2] ) ;
      result[idx3] -= contract_O2O1_2( C1 , C2 , H1H2_degenerate ) ;

      ////////////////// Dimeson -> Dimeson mixing terms
      const size_t idx4 = 36 + loc_idx ;
      
      // O_2 O_2 -- term 1 is positive 
      precompute_block( C1 , bwdH1 , g5[mu1] , L1 , tildeg5[mu2] ) ;
      precompute_block( C2 , bwdH2 , gi[nu1] , L2 , tildegi[nu2] ) ;
      result[idx4]  = contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;

      // O_2 O_2 -- term 2 is -( a b^\dagger )
      precompute_block( C1 , bwdH1 , g5[mu1] , L1 , tildegi[nu2] ) ;
      precompute_block( C2 , bwdH2 , gi[nu1] , L2 , tildeg5[mu2] ) ;  
      result[idx4] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;

      // need to do the others where L1 and L2 are swapped, this is only 
      // a concern for the dimeson - dimeson
      if( L1L2_degenerate == GLU_FALSE ) {
	// O2O2 -- term 3 is -( b a^\dagger )
	precompute_block( C1 , bwdH1 , g5[mu1] , L2 , tildegi[nu2] ) ;
	precompute_block( C2 , bwdH2 , gi[nu1] , L1 , tildeg5[mu2] ) ;  
	result[idx4] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;

	// O2O2 -- term 4 is ( b b^\dagger )
	precompute_block( C1 , bwdH1 , g5[mu1] , L2 , tildeg5[mu2] ) ;
	precompute_block( C2 , bwdH2 , gi[nu1] , L1 , tildegi[nu2] ) ;
	result[idx4] += contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;
      } else {
	result[idx4] *= 2 ;
      }
      //
    }
  }

  free( C1 ) ;
  free( C2 ) ;

  return SUCCESS ;
}
