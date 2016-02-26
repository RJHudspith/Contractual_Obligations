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

// diquark-diquark contractions, now have 4 operators
static int
diquark_diquark( double complex *result ,
		 const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B1 ,
		 const struct spinor B2 ,
		 const struct gamma Cg1 ,
		 const struct gamma Cg2 ,
		 const struct gamma Cg3 ,
		 const struct gamma Cg4 , 
		 const GLU_bool H1H2_degenerate )
{
  // loop of colors
  const size_t Nco = NCNC * NCNC ;

  // counters
  size_t abcd , a , b , c , d ;

  // sum is the two forward propagating diquarks minus the cross term
  register double complex sum1 = 0.0 , sum2 = 0.0 ;

  // temporaries are 
  struct block *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco*sizeof( struct block ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco*sizeof( struct block ) ) != 0 ) {
    goto memfree ;
  }

  // Eq.37 in note is O_1 O_1^\dagger
  {
    sum1 = sum2 = 0.0 ;
    precompute_block( C1 , transpose_spinor( U ) , Cg1 , D , Cg2 ) ;
    precompute_block( C2 , B1 , Cg3 , transpose_spinor( B2 ) , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual forward streaming term
      sum1 += spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) *
	+spinmatrix_trace( C2[ element( c , a , d , b ) ].M ) ;
      // cross color term by interchanging Bs
      sum2 += spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) * 
	spinmatrix_trace( C2[ element( c , b , d , a ) ].M ) ;
    }
    // if the heavies are the same particle we have a cross term
    if( H1H2_degenerate == GLU_TRUE ) {
      result[0] = 4*( sum1 - sum2 ) ;
    } else {
      result[0] = 4*( sum1 ) ;
    }
  }
  // Eq.41 in note is O_1 O_2^\dagger
  {
    sum1 = sum2 = 0.0 ;
    precompute_block( C1 , transpose_spinor( B2 ) , Cg3 , B1 , Cg2 ) ;
    precompute_block( C2 , transpose_spinor( U ) , Cg1 , D , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum1 += 
	trace_prod_spinmatrices( C1[ element( c , a , d , b ) ].M ,
				 C2[ element( a , c , b , d ) ].M ) ;
      sum2 += 
	trace_prod_spinmatrices( C1[ element( c , b , d , a ) ].M ,
				 C2[ element( a , c , b , d ) ].M ) ;
    }
    // if the heavies are the same particle we have a cross term
    if( H1H2_degenerate == GLU_TRUE ) {
      result[1] = 4*( sum1 - sum2 ) ;
    } else {
      result[1] = 4*( sum1 ) ;
    }
  }
  // Eq.43 in note is O_2 O_1^{\dagger}
  {
    sum1 = sum2 = 0.0 ;
    precompute_block( C1 , transpose_spinor( U ) , Cg1 , B1 , Cg4 ) ;
    precompute_block( C2 , transpose_spinor( B2 ) , Cg3 , D , Cg2 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum1 += trace_prod_spinmatrices( C1[ element( a , c , d , b ) ].M ,
				       C2[ element( c , a , b , d ) ].M ) ;
      sum2 += trace_prod_spinmatrices( C1[ element( a , c , d , a ) ].M ,
				       C2[ element( c , b , b , d ) ].M ) ;
    }
    // if the heavies are the same particle we have a cross term
    if( H1H2_degenerate == GLU_TRUE ) {
      result[2] = 4*( sum1 - sum2 ) ;
    } else {
      result[2] = 4*( sum1 ) ;
    }
  }
  // Eq.39 in note is O_2 O_2^\dagger
  {
    sum1 = sum2 = 0.0 ;
    precompute_block( C1 , transpose_spinor( U ) , Cg1 , 
		           transpose_spinor( B2 ) , Cg3 ) ;
    precompute_block( C2 , B1 , Cg2 , D , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum1 += 
	spinmatrix_trace( C1[ element( a , c , d , b ) ].M ) * 
	spinmatrix_trace( C2[ element( c , a , b , d ) ].M ) ;
      sum2 +=
	spinmatrix_trace( C1[ element( a , c , d , a ) ].M ) * 
	spinmatrix_trace( C2[ element( c , b , b , d ) ].M ) ;
    }
    // if the heavies are the same particle we have a cross term
    if( H1H2_degenerate == GLU_TRUE ) {
      result[3] = 4*( sum1 - sum2 ) ;
    } else {
      result[3] = 4*( sum1 ) ;
    }
  }

  free( C1 ) ;
  free( C2 ) ;
  return SUCCESS ;

 memfree :
  free( C1 ) ;
  free( C2 ) ;
  printf( "[TETRA] corr_malloc failure in diquark-diquark\n" ) ;
  return FAILURE ;
}

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
  
  // precompute Cgamma matrices
  const struct gamma Cg5 = CGmu( g5 , GAMMAS ) ;
  const struct gamma Cgi = CGmu( gi , GAMMAS ) ;

  // compute the tilded contributions
  const struct gamma tildeCg5 = gt_Gdag_gt( Cg5 , gt ) ;
  const struct gamma tildeCgi = gt_Gdag_gt( Cgi , gt ) ;

  // poke in the first set of tetra contractions, top left matrix
  if( diquark_diquark( result+0x0 , L1 , L2 , bwdH1 , bwdH2 , 
		       Cg5 , tildeCg5 , Cgi , tildeCgi ,
		       H1H2_degenerate ) == FAILURE ) {
    return FAILURE ;
  }

  // top right matrix is the interchange of 2 -> 4
  if( diquark_diquark( result+0x4 , L1 , L2 , bwdH1 , bwdH2 , 
		       Cg5 , tildeCgi , Cgi , tildeCg5 , 
		       H1H2_degenerate ) == FAILURE ) {
    return FAILURE ;
  }

  // bottom left matrix is the interchange of 1 -> 3
  if( diquark_diquark( result+0x8 , L1 , L2 , bwdH1 , bwdH2 ,
		       Cgi , tildeCg5 , Cg5 , tildeCgi , 
		       H1H2_degenerate ) == FAILURE ) {
    return FAILURE ;
  }

  // here we do the Cg5 <-> Cgi interchanges, bottom right corner of matrix
  if( diquark_diquark( result+0xc , L1 , L2 , bwdH1 , bwdH2 , 
		       Cgi , tildeCgi , Cg5 , tildeCg5 , 
		       H1H2_degenerate ) == FAILURE ) {
    return FAILURE ;
  }

  // compute the meson-meson, if L1 is degenerate with L2 we have exact isospin 
  // and need to flip a sign
  result[ 16 ] = dimeson( L1 , L2 , bwdH1 , bwdH2 , gt , gi , g5 , 
			  L1L2_degenerate , H1H2_degenerate ) ;
  if( L1L2_degenerate == GLU_FALSE ) {
    result[ 16 ] += dimeson( L2 , L1 , bwdH1 , bwdH2 , gt , gi , g5 , 
			     L1L2_degenerate , H1H2_degenerate ) ;
  } else {
    result[ 16 ] *= 2 ;
  }

  // leave if something goes wrong (malloc failure in this case) 
  if( result[16] == sqrt(-1) ) {
    return FAILURE ;
  }

  return SUCCESS ;
}
