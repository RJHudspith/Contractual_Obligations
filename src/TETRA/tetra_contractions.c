/**
   @file tetra_contractions.c
   @brief tetraquark contractions

   TODO :: write some tests for this!
 */
#include "common.h"

#include "contractions.h"       // gamma_mul_r()
#include "gammas.h"             // Cgmu()
#include "spinmatrix_ops.h"     // trace_prod, get_spinmatrix()
#include "spinor_ops.h"         // transpose_spinor()
#include "tetra_contractions.h" // alphabetising

// dimeson contraction
static double complex
dimeson( const struct spinor U , // u prop
	 const struct spinor D , // d prop
	 const struct spinor B1 , // adjoint of B1
	 const struct spinor B2 , // adjoint of B2
	 const struct gamma gt ,
	 const struct gamma gi ,
	 const struct gamma g5 , 
	 const GLU_bool L1L2_degenerate ,
	 const GLU_bool H1H2_degenerate )
{
  // flattened number of colors we loop over
  const size_t Nco = NCNC * NCNC ;

  // index loops
  size_t abcd , a , b , c , d ;

  // temporaries blocks of spinmatrix data
  struct block *C1 = NULL , *C2 = NULL ;

  // precompute our tilded temporaries
  const struct gamma tildeg5 = gt_Gdag_gt( g5 , gt )  ;
  const struct gamma tildegi = gt_Gdag_gt( gi , gt )  ;

  // sums
  register double complex sum1 = 0.0 ;
  register double complex sum2 = 0.0 ;
  register double complex sum3 = 0.0 ;
  register double complex sum4 = 0.0 ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco * sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco * sizeof( struct block ) ) != 0 ) {
    goto memfree ;
  }

  // first term is (O_3 O_3^\dagger)^{(1)} = Eq.46 in note
  {
    precompute_block( C1 , B1 , g5 , U , tildeg5 ) ;
    precompute_block( C2 , B2 , gi , D , tildegi ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum1 += 
	spinmatrix_trace( C1[ element( c , a , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , b , d ) ].M ) ;

      // cross term
      sum3 += 
	spinmatrix_trace( C1[ element( c , b , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , b , d ) ].M ) ;
    }
  }
  // second part is (O_3 O_3^\dagger)^{(2)} Eq.48 in note
  {
    precompute_block( C1 , B1 , gi , D , tildeg5 ) ;
    precompute_block( C2 , B2 , g5 , U , tildegi ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum2 += 
	spinmatrix_trace( C1[ element( c , b , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , a , d ) ].M ) ;
      // cross term
      sum4 += 
	spinmatrix_trace( C1[ element( c , a , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , a , d ) ].M ) ;
    }
  }

  free( C1 ) ;
  free( C2 ) ;

  // if the heavies are degenerate we include the cross terms
  if( H1H2_degenerate == GLU_TRUE ) {
    return ( sum1 - sum2 ) - ( sum3 - sum4 ) ;
  } else {
    return ( sum1 - sum2 ) ;
  }

 memfree :
  free( C1 ) ;
  free( C2 ) ;
  printf( "[TETRA] corr_malloc failure in dimeson\n" ) ;
  return sqrt(-1) ;
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
  if( H1H2_degenerate == GLU_TRUE ) {
    return 4*( sum1 - sum2 ) ;
  } else {
    return 4*( sum1 ) ;
  }
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

  // loop of colors
  const size_t Nco = NCNC * NCNC ;

  // temporaries are 
  struct block *C1 = NULL , *C2 = NULL , *C3 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco*sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco*sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C3 , 16 , Nco*sizeof( struct block ) ) != 0 ) {
  }

  // precompute Cgamma matrices and Pull out the heavy contribution
  // which is the block of bs with a gamma_i inside
  const struct gamma Cgi = CGmu( gi , GAMMAS ) ;
  const struct gamma tildeCgi = gt_Gdag_gt( Cgi , gt ) ;
  precompute_block( C2 , bwdH1 , Cgi , transpose_spinor( bwdH2 ) , tildeCgi ) ;

  // Operator
  //
  //  ( u^T C \gamma_5 d ) ( \bar{b} C \gamma_i \bar{b}^T ) 
  //
  const struct gamma Cg5 = CGmu( g5 , GAMMAS ) ;
  const struct gamma tildeCg5 = gt_Gdag_gt( Cg5 , gt ) ;
  precompute_block( C1 , transpose_spinor( L1 ) , Cg5 , L2 , tildeCg5 ) ;
  result[0] = contract_O1O1( C1 , C2 , H1H2_degenerate ) ;

  // Contraction
  //
  //  ( u^T C d )
  //
  const struct gamma C = CGmu( GAMMAS[ IDENTITY ] , GAMMAS ) ;
  const struct gamma tildeC = gt_Gdag_gt( C , gt ) ;
  precompute_block( C1 , transpose_spinor( L1 ) , C , L2 , tildeC ) ;

  // Contraction
  //
  //  ( u^T C \gamma_i \gamma_5 d )
  //
  const struct gamma Cgig5 = CGmu( GAMMAS[ GAMMA_5 + mu + 1 ]  , GAMMAS ) ;
  const struct gamma tildeCgig5 = gt_Gdag_gt( Cgig5 , gt ) ;
  precompute_block( C3 , transpose_spinor( L1 ) , Cgig5 , L2 , tildeCgig5 ) ;

  // and the weird ones
  size_t j , count = 0 ;
  for( j = 0 ; j < ND-1 ; j++ ) { 
    if( j == mu ) continue ;

    // precompute the Bs first as the result can be used twice
    struct gamma gigj ;
    gamma_mmul( &gigj , gi , GAMMAS[ j ] ) ;
    const struct gamma Cgigj = CGmu( gigj , GAMMAS ) ;
    const struct gamma tildeCgigj = gt_Gdag_gt( Cgigj , gt ) ;
    precompute_block( C2 , bwdH1 , Cgigj , 
		      transpose_spinor( bwdH2 ) , tildeCgigj ) ;

    result[ 1 + count ] = contract_O1O1( C1 , C2 , H1H2_degenerate ) ;
    result[ 2 + count ] = contract_O1O1( C3 , C2 , H1H2_degenerate ) ;

    count+=2 ;
  }

  // free these temporaries
  free( C1 ) ;
  free( C2 ) ;
  free( C3 ) ;

  // dimeson contractions
  result[ TETRA_NOPS - 1 ] = dimeson( L1 , L2 , bwdH1 , bwdH2 , gt , gi , g5 , 
				      L1L2_degenerate , H1H2_degenerate ) ;
  if( L1L2_degenerate == GLU_FALSE ) {
    result[ TETRA_NOPS - 1 ] += dimeson( L2 , L1 , bwdH1 , bwdH2 , gt , gi , g5 , 
					 L1L2_degenerate , H1H2_degenerate ) ;
  } else {
    result[ TETRA_NOPS - 1 ] *= 2 ;
  }

  // leave if something goes wrong (malloc failure in this case) 
  if( result[ TETRA_NOPS - 1 ] == sqrt(-1) ) {
    return FAILURE ;
  }

  return SUCCESS ;
}
