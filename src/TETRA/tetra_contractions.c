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
int
diquark_diquark( double complex result[ TETRA_NOPS-1 ] ,
		 const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B , // full adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu )
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

  const struct gamma Cg5      = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  const struct gamma tildeCg5 = gt_Gdag_gt( Cg5 , GAMMAS )  ;
  const struct gamma Cgi      = CGmu( GAMMAS[ mu ] , GAMMAS ) ;
  const struct gamma tildeCgi = gt_Gdag_gt( Cgi , GAMMAS ) ;

  // Eq.37 in note is O_1 O_1^\dagger
  {
    sum = 0.0 ;
    precompute_block( C1 , U , Cg5 , D , tildeCg5 ) ;
    precompute_block( C2 , B , Cgi , B , tildeCgi ) ;
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
    precompute_block( C1 , B , Cgi , B , tildeCg5 ) ;
    precompute_block( C2 , U , Cg5 , D , tildeCgi ) ;
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
    precompute_block( C1 , U , Cg5 , B , tildeCgi ) ;
    precompute_block( C2 , B , Cgi , D , tildeCg5 ) ;
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
    precompute_block( C1 , U , Cg5 , B , tildeCg5 ) ;
    precompute_block( C2 , B , Cgi , D , tildeCgi ) ;
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
double complex
dimeson( const struct spinor U , // u prop
	 const struct spinor D , // d prop
	 const struct spinor B , // adjoint of B
	 const struct gamma *GAMMAS ,
	 const size_t mu )
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

  // precompute our temporaries
  const struct gamma g5      = GAMMAS[ GAMMA_5 ] ;
  const struct gamma tildeg5 = gt_Gdag_gt( g5 , GAMMAS )  ;
  const struct gamma gi      = GAMMAS[ mu ] ;
  const struct gamma tildegi = gt_Gdag_gt( gi , GAMMAS )  ;

  register double complex sum = 0.0 ;
  // first term is (O_3 O_3^\dagger)^{(1)} = Eq.46 in note
  {
    precompute_block( C1 , B , g5 , U , tildegi ) ;
    precompute_block( C2 , B , gi , D , tildeg5 ) ;
    for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum += 
	spinmatrix_trace( C1[ element( c , a , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , b , d ) ].M ) ;

      // cross term
      sum -= 
	spinmatrix_trace( C1[ element( c , b , a , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , b , d ) ].M ) ;
    }
  }
  // second part is (O_3 O_3^\dagger)^{(2)} Eq.48 in note, global minus on this one!
  {
    precompute_block( C1 , B , gi , D , tildeg5 ) ;
    precompute_block( C2 , B , g5 , U , tildegi ) ;
    for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
      get_abcd( &a , &b , &c , &d , abcd ) ;
      // usual meson product
      sum -= 
	spinmatrix_trace( C1[ element( c , b , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , a , a , d ) ].M ) ;
      // cross term
      sum += 
	spinmatrix_trace( C1[ element( c , a , b , c ) ].M ) *
	spinmatrix_trace( C2[ element( d , b , a , d ) ].M ) ;
    }
  }

  free( C1 ) ;
  free( C2 ) ;
  return sum ;

 memfree :
  free( C1 ) ;
  free( C2 ) ;
  printf( "[TETRA] corr_malloc failure in dimeson\n" ) ;
  return sqrt(-1) ;
}
