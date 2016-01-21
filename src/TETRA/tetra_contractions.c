/**
   @file tetra_contractions.c
   @brief tetraquark contractions

   TODO :: write some tests for this!
 */
#include "common.h"

#include "contractions.h"   // gamma_mul_r()
#include "gammas.h"         // Cgmu()
#include "spinmatrix_ops.h" // trace_prod_spinmatrices()

// little storage struct
struct cc {
  double complex M[ NSNS ] ;
} ;

// this is a really out of order memory access, should be in spinor_ops no?
static void
get_spinmatrix( double complex *s , 
		const struct spinor S ,
		const size_t c1 , 
		const size_t c2 )
{
  // pokes specific colors c1,c2 out from the spinor
  // each spinor is 4x4x3x3 chooses a component of the 3x3
  // color matrix and hence gives a 4x4 spinmatrix
  const double complex *dS = (const double complex*)( S.D ) ; 
  dS += c2 + NC * c1  ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    *s = *dS ; s++ ; dS += NCNC ;
  }
  return ;
}

// get an element from a b c d
static size_t
element( const size_t a , const size_t b , const size_t c , const size_t d ) 
{
  return d + NC * ( c + NC * ( b + NC * a ) ) ;
}

// little function to determine a,b,c and d. d runs fastest, a slowest
static void
get_abcd( size_t *a , size_t *b , size_t *c , size_t *d , const size_t abcd )
{
  // cache this
  const size_t NCNCNC = NCNC * NC ;
  *a = ( abcd ) / ( NCNCNC ) ;
  *b = ( abcd % NCNCNC ) / ( NCNC ) ;
  *c = ( abcd % NCNC ) / ( NC ) ;
  *d = ( abcd % NC ) ;
}

// precompute a "block" is a spinmatrix
static void
precompute_block( struct cc *C1 ,
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

// diquark-diquark contractions, now have 4 operators
static int
diquark_diquark( double complex result[ TETRA_NOPS-1 ] ,
		 const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B ,
		 const struct gamma Cg1 ,
		 const struct gamma Cg2 ,
		 const struct gamma Cg3 ,
		 const struct gamma Cg4 )
{
  // loop of colors
  const size_t Nco = NCNC * NCNC ;

  // counters
  size_t abcd , a , b , c , d ;

  // sum is the two forward propagating diquarks minus the cross term
  register double complex sum = 0.0 ;

  // temporaries are 
  struct cc *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco*sizeof( struct cc ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco*sizeof( struct cc ) ) != 0 ) {
    goto memfree ;
  }

  // Eq.37 in note is O_1 O_1^\dagger
  {
    sum = 0.0 ;
    precompute_block( C1 , U , Cg1 , D , Cg2 ) ;
    precompute_block( C2 , B , Cg3 , B , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum += 
	spinmatrix_trace( C1[ element( a , c , b , d ) ].M ) * 
	( spinmatrix_trace( C2[ element( c , a , d , b ) ].M ) - 
	  spinmatrix_trace( C2[ element( c , b , d , a ) ].M ) ) ;
    }
    result[0] = sum ;
  }
  // Eq.41 in note is O_1 O_2^\dagger
  {
    sum = 0.0 ;
    precompute_block( C1 , B , Cg3 , B , Cg2 ) ;
    precompute_block( C2 , U , Cg1 , D , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum += 
	trace_prod_spinmatrices( C1[ element( c , a , d , b ) ].M ,
				 C2[ element( a , c , b , d ) ].M ) -
	trace_prod_spinmatrices( C1[ element( c , b , d , a ) ].M ,
				 C2[ element( a , c , b , d ) ].M ) ;
    }
    result[1] = sum ;
  }
  // Eq.43 in note is O_2 O_1^{\dagger}
  {
    sum = 0.0 ;
    precompute_block( C1 , U , Cg1 , B , Cg4 ) ;
    precompute_block( C2 , B , Cg3 , D , Cg2 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum += 
	trace_prod_spinmatrices( C1[ element( a , c , d , b ) ].M ,
				 C2[ element( c , a , b , d ) ].M ) -
	trace_prod_spinmatrices( C1[ element( a , c , d , a ) ].M ,
				 C2[ element( c , b , b , d ) ].M ) ;
    }
    result[2] = sum ;
  }
  // Eq.39 in note is O_2 O_2^\dagger
  {
    sum = 0.0 ;
    precompute_block( C1 , U , Cg1 , B , Cg2 ) ;
    precompute_block( C2 , B , Cg3 , D , Cg4 ) ;
    for( abcd = 0 ; abcd < Nco ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      sum += 
	// block1
	spinmatrix_trace( C1[ element( a , c , d , b ) ].M ) * 
	spinmatrix_trace( C2[ element( c , a , b , d ) ].M ) - 
	// block2
	spinmatrix_trace( C1[ element( a , c , d , a ) ].M ) * 
	spinmatrix_trace( C2[ element( c , b , b , d ) ].M ) ;
    }
    result[3] = sum ;
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
	 const struct spinor B , // adjoint of B
	 const struct gamma gt ,
	 const struct gamma gi ,
	 const struct gamma g5 , 
	 const GLU_bool L1L2_degenerate )
{
  // flattened number of colors we loop over
  const size_t Nco = NCNC * NCNC ;

  // index loops
  size_t abcd , a , b , c , d ;

  // temporaries blocks of spinmatrix data
  struct cc *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco * sizeof( struct cc ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco * sizeof( struct cc ) ) != 0 ) {
    goto memfree ;
  }

  // precompute our tilded temporaries
  const struct gamma tildeg5 = gt_Gdag_gt( g5 , gt )  ;
  const struct gamma tildegi = gt_Gdag_gt( gi , gt )  ;

  register double complex sum1 = 0.0 ;
  // first term is (O_3 O_3^\dagger)^{(1)} = Eq.46 in note
  {
    precompute_block( C1 , B , g5 , U , tildegi ) ;
    precompute_block( C2 , B , gi , D , tildeg5 ) ;
    for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum1 += 
	spinmatrix_trace( C1[ element( c , a , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , b , d ) ].M ) ;

      // cross term
      sum1 -= 
	spinmatrix_trace( C1[ element( c , b , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , b , d ) ].M ) ;
    }
  }
  register double complex sum2 = 0.0 ;
  // second part is (O_3 O_3^\dagger)^{(2)} Eq.48 in note
  {
    precompute_block( C1 , B , gi , D , tildeg5 ) ;
    precompute_block( C2 , B , g5 , U , tildegi ) ;
    for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum2 += 
	spinmatrix_trace( C1[ element( c , b , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , a , d ) ].M ) ;
      // cross term
      sum2 -= 
	spinmatrix_trace( C1[ element( c , a , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , a , d ) ].M ) ;
    }
  }

  free( C1 ) ;
  free( C2 ) ;

  // if u and d are the same, add them
  if( L1L2_degenerate == GLU_TRUE ) {
    return sum1 + sum2 ;
  } else {
    return sum1 - sum2 ;
  }

 memfree :
  free( C1 ) ;
  free( C2 ) ;
  printf( "[TETRA] corr_malloc failure in dimeson\n" ) ;
  return sqrt(-1) ;
}

// perform the contraction
int
tetras( double complex result[ TETRA_NOPS ] ,
	const struct spinor L1 , 
	const struct spinor L2 ,
	const struct spinor bwdH ,
	const struct gamma *GAMMAS ,
	const size_t mu , 
	const GLU_bool L1L2_degenerate )
{
  // rename gamma matrices to match the contractions above
  const struct gamma gt = GAMMAS[ GAMMA_3 ] ;
  const struct gamma gi = GAMMAS[ mu ] ;
  const struct gamma g5 = GAMMAS[ GAMMA_5 ] ;
  
  // precompute Cgamma matrices
  const struct gamma Cg5 = CGmu( g5 , GAMMAS ) ;
  const struct gamma Cgi = CGmu( gi , GAMMAS ) ;

  // compute the tilded contributions
  const struct gamma tildeCg5 = gt_Gdag_gt( Cg5 , gt ) ;
  const struct gamma tildeCgi = gt_Gdag_gt( Cgi , gt ) ;

  // poke in the first set of tetra contractions, top left matrix
  if( diquark_diquark( result , L1 , L2 , bwdH , 
		       Cg5 , tildeCg5 , Cgi , tildeCgi ) == FAILURE ) {
    return FAILURE ;
  }

  // top right matrix is the interchange of 2 -> 4
  if( diquark_diquark( result+4 , L1 , L2 , bwdH , 
		       Cg5 , tildeCgi , Cgi , tildeCg5 ) == FAILURE ) {
    return FAILURE ;
  }

  // bottom left matrix is the interchange of 1 -> 3
  if( diquark_diquark( result+8 , L1 , L2 , bwdH , 
		       Cgi , tildeCg5 , Cg5 , tildeCgi ) == FAILURE ) {
    return FAILURE ;
  }

  // here we do the Cg5 <-> Cgi interchanges, bottom right corner of matrix
  if( diquark_diquark( result+12 , L1 , L2 , bwdH , 
		       Cgi , tildeCgi , Cg5 , tildeCg5 ) == FAILURE ) {
    return FAILURE ;
  }

  // compute the meson-meson, if L1 is degenerate with L2 we have exact isospin 
  // and need to flip a sign
  result[ 16 ] = dimeson( L1 , L2 , bwdH , gt , gi , g5 , L1L2_degenerate ) ;
  if( L1L2_degenerate == GLU_FALSE ) {
    result[ 16 ] += dimeson( L2 , L1 , bwdH , gt , gi , g5 , L1L2_degenerate ) ;
  } else {
    result[ 16 ] *= 2 ;
  }
  // leave if something goes wrong (malloc failure in this case) 
  if( result[16] == sqrt(-1) ) {
    return FAILURE ;
  }

  return SUCCESS ;
}
