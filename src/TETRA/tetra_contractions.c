/**
   @file tetra_contractions.c
   @brief tetraquark contractions

   TODO: we should speed this up! can precompute all gamma
   multiplies. Could break up NC loops and do precomputations.
   Can we do all the multiplies outside of the abcd loop?
 */
#include "common.h"

#include "contractions.h"   // gamma_mul_r()
#include "gammas.h"         // Cgmu()
#include "spinmatrix_ops.h" // trace_prod_spinmatrices()

// little storage struct
struct cc {
  double complex M[ NSNS ] ;
} ;

// this is a really out of order memory access
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

// little function to determine a,b,c and d
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

// precomputes some spinmatrices
static void
precompute_C1C2( struct cc *C1 ,
		 struct cc *C2 ,
		 const struct spinor S1 ,
		 const struct spinor S2 ,
		 const struct spinor S3 ,
		 const struct spinor S4 ,
		 const struct gamma G12 ,
		 const struct gamma G34 ,
		 const struct gamma *GAMMAS )
{
  struct spinor S1g12 = S1 , S2g12dag = S2 ;
  const struct gamma G12dag = gt_Gdag_gt( G12 , GAMMAS ) ;
  gamma_mul_r( &S1g12    , G12    ) ;
  gamma_mul_r( &S2g12dag , G12dag ) ;

  struct spinor S3g34 = S3 , S4g34dag = S4 ;
  const struct gamma G34dag = gt_Gdag_gt( G34 , GAMMAS ) ;
  gamma_mul_r( &S3g34    , G34    ) ;
  gamma_mul_r( &S4g34dag , G34dag ) ;

  double complex D1[ NSNS ] , D2[ NSNS ] , D3[ NSNS ] , D4[ NSNS ] ;
  size_t ab , cd ;
  for( ab = 0 ; ab < NCNC ; ab++ ) {

    const size_t a = ab / NC ;
    const size_t b = ab % NC ;

    get_spinmatrix( D1 , S1g12    , a , b ) ;
    get_spinmatrix( D3 , S3g34    , a , b ) ;

    for( cd = 0 ; cd < NCNC ; cd++ ) {

      const size_t c = cd / NC ;
      const size_t d = cd % NC ;

      get_spinmatrix( D2 , S2g12dag , c , d ) ;
      get_spinmatrix( D4 , S4g34dag , c , d ) ;

      spinmatrix_multiply( C1 -> M , D1 , D2 ) ;
      spinmatrix_multiply( C2 -> M , D3 , D4 ) ;

      C1++ , C2++ ;
    }
  }

  return ;
}

// diquark-diquark contraction
double complex
diquark_diquark( const struct spinor U ,
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

  // temporaries are UCg5.DCg5dag and BCgmu.BCgmudag
  precompute_C1C2( C1 , C2 , 
		   U , D , B , B ,
		   CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ,
		   CGmu( GAMMAS[ mu ] , GAMMAS ) , 
		   GAMMAS ) ;

  // loop color indices again
  for( abcd = 0 ; abcd < Nco ; abcd++ ) {
    // get the color indices from linearised index abcd
    get_abcd( &a , &b , &c , &d , abcd ) ;

    // compute the sum
    sum += 
      // TrS[ U_{ac}.Cg5.D_{bd}.Cg5T.B_{ac}.Cgi.B_{bd}.CgiT ]
      +trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M , 
				C2[ element( a , c , b , d ) ].M )
      // TrS[ U_{ac}.Cg5.D_{bd}.Cg5T.B_{ad}.Cgi.B_{bc}.CgiT ]
      -trace_prod_spinmatrices( C1[ element( a , c , b , d ) ].M , 
			       C2[ element( a , d , b , c ) ].M ) ;
  }

  free( C1 ) ;
  free( C2 ) ;
  return sum ;

 memfree :
  free( C1 ) ;
  free( C2 ) ;
  printf( "[TETRA] corr_malloc failure in diquark-diquark\n" ) ;
  return sqrt(-1) ;
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

  // temporaries are Bg5 . Ug5dag and Bgmu . Dgmudag
  struct cc *C1 = NULL , *C2 = NULL ;

  // make everything a NaN if we fail to allocate memory
  if( corr_malloc( (void**)&C1 , 16 , Nco * sizeof( struct cc ) ) != 0 || 
      corr_malloc( (void**)&C2 , 16 , Nco * sizeof( struct cc ) ) != 0 ) {
    goto memfree ;
  }

  // precompute our temporaries
  precompute_C1C2( C1 , C2 , 
		   B , U , B , D ,
		   GAMMAS[ GAMMA_5 ] ,
		   GAMMAS[ mu ] , 
		   GAMMAS ) ;

  register double complex sum = 0.0 ;
  // loop colors
  for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
    // get the color indices from linearised index abcd
    get_abcd( &a , &b , &c , &d , abcd ) ;

    // sum the four possible terms
    sum += 
      // +TrS[ B_{ac}.g5.U_{ac}.g5T.B_{bd}.gi.B_{bd}.giT ]
      +trace_prod_spinmatrices( C1[ element( a , c , a , c ) ].M , 
				C2[ element( b , d , b , d ) ].M ) 
      // -TrS[ B_{ad}.g5.U_{ac}.g5T.B_{bc}.gi.B_{bd}.giT ]
      -trace_prod_spinmatrices( C1[ element( a , d , a , c ) ].M , 
				C2[ element( b , c , b , d ) ].M ) 
      // -TrS[ B_{ac}.g5.U_{ad}.g5T.B_{bd}.gi.B_{bc}.giT ]
      -trace_prod_spinmatrices( C1[ element( a , c , a , d ) ].M , 
				C2[ element( b , d , b , c ) ].M ) 
      // +TrS[ B_{ad}.g5.U_{ad}.g5T.B_{bc}.gi.B_{bc}.giT ]
      +trace_prod_spinmatrices( C1[ element( a , d , a , d ) ].M , 
				C2[ element( b , c , b , c ) ].M ) ;
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
