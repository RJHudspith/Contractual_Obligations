/**
   @file tetra_contractions.c
   @brief tetraquark contractions

   TODO :: perhaps move some of these to LINALG?
 */
#include "common.h"

#include "gammas.h"         // Cgmu()
#include "spinmatrix_ops.h" // trace_prod_spinmatrices()

// this is a really out of order memory access
static void
get_spinmatrix( double complex *s , 
		const struct spinor S ,
		const size_t c1 , 
		const size_t c2 )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      s[ d2 + NS * d1 ] = S.D[ d1 ][ d2 ].C[ c1 ][ c2 ] ;
    }
  }
  return ;
}

// this is really out of order, in general
static void
get_spinmatrix_mulgamma( double complex *s , 
			 const struct spinor S ,
			 const size_t c1 , 
			 const size_t c2 ,
			 const struct gamma G )
{
  get_spinmatrix( s , S , c1 , c2 ) ;
  spinmatrix_gamma( s , G ) ;
  return ;
}

// read as Tr[ (s1_{c1,c2}.G1).(s2_{c3,c4}.G2).(s3_{c5,c6}.G3).(s4_{c7,c8}.G4) ]
// i.e. Trace left to right product of spinmatrices with color components c_i
static double complex
c4( const struct spinor s1 , const size_t c1 , const size_t c2 , const struct gamma G1 ,
    const struct spinor s2 , const size_t c3 , const size_t c4 , const struct gamma G2 ,
    const struct spinor s3 , const size_t c5 , const size_t c6 , const struct gamma G3 ,
    const struct spinor s4 , const size_t c7 , const size_t c8 , const struct gamma G4 )
{
  double complex D1[ NSNS ] , D2[ NSNS ] , prod1[ NSNS ] , prod2[ NSNS ] ;
  get_spinmatrix_mulgamma( D1 , s1 , c1 , c2 , G1 ) ; 
  get_spinmatrix_mulgamma( D2 , s2 , c3 , c4 , G2 ) ; 
  spinmatrix_multiply( prod1 , D1 , D2 ) ;
  get_spinmatrix_mulgamma( D1 , s3 , c5 , c6 , G3 ) ; 
  get_spinmatrix_mulgamma( D2 , s4 , c7 , c8 , G4 ) ; 
  spinmatrix_multiply( prod2 , D1 , D2 ) ;
  return trace_prod_spinmatrices( prod1 , prod2 ) ;
}

// diquark-diquark contraction
double complex
diquark_diquark( const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B , // full adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu )
{
  // precompute Cg5 and Cg5T
  const struct gamma Cg5  = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  const struct gamma Cg5T = CGmuT( Cg5 , GAMMAS ) ;

  // ang Cgmu
  const struct gamma Cgmu = CGmu( GAMMAS[ mu ] , GAMMAS ) ;
  const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;

  size_t a , b ;
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  for( a = 0 ; a < NC ; a++ ) {
    for( b = 0 ; b < NC ; b++ ) {
      // [ Uaa.Cg5.Dbb.Cg5T.Baa.Cgmu.Bbb.CgmuT ]
      sum1 += c4( U , a , a , Cg5 ,  
		  D , b , b , Cg5T ,
		  B , a , a , Cgmu , 
		  B , b , b , CgmuT ) ;
      // [ Uaa.Cg5.Dbb.Cg5T.Baa.Cgmu.Bbb.CgmuT ]
      sum2 += c4( U , a , a , Cg5 ,  
		  D , b , b , Cg5T ,
		  B , a , b , Cgmu , 
		  B , b , a , CgmuT ) ;
    }
  }
  return sum1 - sum2 ;
}

// diquark-diquark contraction
double complex
dimeson_dimeson( const struct spinor U ,  // u prop
		 const struct spinor D ,  // d prop
		 const struct spinor B , // adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu )
{
  // four possible sums
  register double complex sum1 = 0.0 , sum2 = 0.0 , 
    sum3 = 0.0 , sum4 = 0.0 ;

  // cache the common gamma matrices
  const struct gamma g5 = GAMMAS[ GAMMA_5 ] , gi = GAMMAS[ mu ] ;

  size_t a , b ;
  for( a = 0 ; a < NC ; a++ ) {
    for( b = 0 ; b < NC ; b++ ) {
      // [ Baa.g5.Uaa.gi.Bbb.g5.Dbb.gi ]
      sum1 += c4( B , a , a , g5 , 
		  U , a , a , gi ,
		  B , b , b , g5 , 
		  D , b , b , gi ) ;
      // [ Bab.g5.Uaa.gi.Bba.g5.Dbb.gi ]
      sum2 += c4( B , a , b , g5 , 
		  U , a , a , gi ,
		  B , b , a , g5 , 
		  D , b , b , gi ) ;
      // [ Baa.g5.Uab.gi.Bbb.g5.Dba.gi ]
      sum3 += c4( B , a , a , g5 , 
		  U , a , b , gi ,
		  B , b , b , g5 , 
		  D , b , a , gi ) ;
      // [ Bab.g5.Uab.gi.Bba.g5.Dba.gi ]
      sum4 += c4( B , a , b , g5 , 
		  U , a , b , gi ,
		  B , b , a , g5 , 
		  D , b , a , gi ) ;
    }
  }
  return sum1 - sum2 - sum3 + sum4 ;
}
