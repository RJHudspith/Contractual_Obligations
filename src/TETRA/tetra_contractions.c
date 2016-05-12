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
    sum1 += spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) *
      spinmatrix_trace( C2[ element( c , a , d , b ) ].M ) ;
    // cross color term by interchanging Bs
    sum2 += spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) * 
      spinmatrix_trace( C2[ element( c , b , d , a ) ].M ) ;
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
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , a , d , b ) ].M ) ;
    // cross term
    sum2 +=
      trace_prod_spinmatrices( C1[ element( a , d , b , c ) ].M ,
			       C2[ element( c , b , d , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*( sum1 - sum2 ) : 4*sum1 ;
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
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , a , c , b ) ].M ) ;
    // cross term
    sum2 +=
      trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M ,
			       C2[ element( d , b , c , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*( sum1 - sum2 ) : 4*sum1 ;
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
      trace_prod_spinmatrices( C1[ element( c , a , b , d ) ].M ,
			       C2[ element( a , c , d , b ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*( sum1 - sum2 ) : 4*sum1 ;
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
      trace_prod_spinmatrices( C1[ element( c , b , a , d ) ].M ,
			       C2[ element( b , c , d , a ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? 4*( sum1 - sum2 ) : 4*sum1 ;
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
      spinmatrix_trace( C1[ element( c , b , a , c ) ].M ) *
      spinmatrix_trace( C2[ element( d , a , b , d ) ].M ) ;
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
      spinmatrix_trace( C1[ element( c , b , b , c ) ].M ) *
      spinmatrix_trace( C2[ element( d , a , a , d ) ].M ) ;
    // cross term
    sum2 += 
      spinmatrix_trace( C1[ element( c , a , b , c ) ].M ) *
      spinmatrix_trace( C2[ element( d , b , a , d ) ].M ) ;
  }
  // if the heavies are the same particle we have a cross term
  return ( H1H2_degenerate == GLU_TRUE ) ? ( sum1 - sum2 ) : sum1 ;
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
  // rename gamma matrices to match the contractions above
  const struct gamma gt = GAMMAS[ GAMMA_T ] ;
  const struct gamma gi = GAMMAS[ mu ] ;
  const struct gamma g5 = GAMMAS[ GAMMA_5 ] ;
  const struct gamma tildegi = gt_Gdag_gt( gi , gt ) ;
  const struct gamma tildeg5 = gt_Gdag_gt( g5 , gt ) ;

  // charge conjugation matrices
  const struct gamma Cgi = CGmu( gi , GAMMAS ) ;
  const struct gamma Cg5 = CGmu( g5 , GAMMAS ) ;
  const struct gamma tildeCgi = gt_Gdag_gt( Cgi , gt ) ;
  const struct gamma tildeCg5 = gt_Gdag_gt( Cg5 , gt ) ;

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

  // O_1 O_1
  precompute_block( C1 , L1T , Cg5 , L2 , tildeCg5 ) ;
  precompute_block( C2 , bwdH1 , Cgi , bwdH2T , tildeCgi ) ;
  result[0] = contract_O1O1( C1 , C2 , H1H2_degenerate ) ;

  // O_1 O_2
  precompute_block( C1 , L1T , Cg5 , L2 , gamma_transpose( tildegi ) ) ;
  precompute_block( C2 , bwdH1 , Cgi , bwdH2T , tildeg5 ) ;
  result[1]  = contract_O1O2_1( C1 , C2 , H1H2_degenerate ) ;

  precompute_block( C1 , L1T , Cg5 , L2 , gamma_transpose( tildeg5 ) ) ;
  precompute_block( C2 , bwdH1 , Cgi , bwdH2T , tildegi ) ;
  result[1] -= contract_O1O2_2( C1 , C2 , H1H2_degenerate ) ;

  // O_2 O_1
  precompute_block( C1 , bwdH1 , gamma_transpose( gi ) , L2 , tildeCg5 ) ;
  precompute_block( C2 , L1T , g5 , bwdH2T , tildeCgi ) ;
  result[2]  = contract_O2O1_1( C1 , C2 , H1H2_degenerate ) ;

  precompute_block( C1 , bwdH1 , gamma_transpose( g5 ) , L2 , tildeCg5 ) ;
  precompute_block( C2 , L1T , gi , bwdH2T , tildeCgi ) ;
  result[2] -= contract_O2O1_2( C1 , C2 , H1H2_degenerate ) ;

  // O_2 O_2
  precompute_block( C1 , bwdH1 , g5 , L1 , tildeg5 ) ;
  precompute_block( C2 , bwdH2 , gi , L2 , tildegi ) ;
  result[3]  = contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;

  precompute_block( C1 , bwdH1 , gi , L2 , tildeg5 ) ;
  precompute_block( C2 , bwdH2 , g5 , L1 , tildegi ) ;  
  result[3] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;

  // need to do the others where L1 and L2 are swapped
  if( L1L2_degenerate == GLU_FALSE ) {
    struct spinor L2T = transpose_spinor( L2 ) ;
    // O1O2
    precompute_block( C1 , L2T , Cg5 , L1 , gamma_transpose( tildegi ) ) ;
    precompute_block( C2 , bwdH1 , Cgi , bwdH2T , tildeg5 ) ;
    result[1] += contract_O1O2_1( C1 , C2 , H1H2_degenerate ) ;

    precompute_block( C1 , L2T, Cg5 , L1 , tildeg5 ) ;
    precompute_block( C2 , bwdH1 , Cgi , bwdH2T , tildegi ) ;
    result[1] -= contract_O1O2_2( C1 , C2 , H1H2_degenerate ) ;

    // O2O1
    precompute_block( C1 , bwdH1 , gamma_transpose( gi ) , L1 , tildeCg5 ) ;
    precompute_block( C2 , L2T , g5 , bwdH2T , tildeCgi ) ;
    result[2] += contract_O2O1_1( C1 , C2 , H1H2_degenerate ) ;

    precompute_block( C1 , bwdH1 , gamma_transpose( g5 ) , L1 , tildeCg5 ) ;
    precompute_block( C2 , L2T , gi , bwdH2T , tildeCgi ) ;
    result[2] -= contract_O2O1_2( C1 , C2 , H1H2_degenerate ) ;

    // O2O2
    precompute_block( C1 , bwdH1 , g5 , L2 , tildeg5 ) ;
    precompute_block( C2 , bwdH2 , gi , L1 , tildegi ) ;
    result[3] += contract_O2O2_1( C1 , C2 , H1H2_degenerate ) ;

    precompute_block( C1 , bwdH1 , gi , L1 , tildeg5 ) ;
    precompute_block( C2 , bwdH2 , g5 , L2 , tildegi ) ;  
    result[3] -= contract_O2O2_2( C1 , C2 , H1H2_degenerate ) ;
  } else {
    result[1] *= 2 ; result[2] *= 2 ; result[3] *= 2 ;
  }

  free( C1 ) ;
  free( C2 ) ;

  return SUCCESS ;
}
