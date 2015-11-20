/**
   @file tetra_contractions.c
   @brief tetraquark contractions
 */
#include "common.h"

#include "contractions.h"   // gamma_mul_r()
#include "gammas.h"         // Cgmu()
#include "spinmatrix_ops.h" // trace_prod_spinmatrices()

// this is a really out of order memory access
static void
get_spinmatrix( double complex *s , 
		const struct spinor S ,
		const size_t c1 , 
		const size_t c2 )
{
  const double complex *dS = (const double complex*)S.D ;
  dS += c2 + NC * c1 ;
  size_t d1d2 ;
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    *s = *dS ; s++ ; dS += NCNC ;
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

// diquark-diquark contraction
double complex
diquark_diquark( const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B , // full adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu )
{
  // precompute Cg5 and Cg5T
  const struct gamma Cg5   = CGmu( GAMMAS[ GAMMA_5 ] , GAMMAS ) ;
  const struct gamma Cg5T  = gt_Gdag_gt( Cg5 , GAMMAS ) ;

  // ang Cgmu
  const struct gamma Cgmu  = CGmu( GAMMAS[ mu ] , GAMMAS ) ;
  const struct gamma CgmuT = gt_Gdag_gt( Cgmu , GAMMAS ) ;

  // loop all colors
  size_t abcd , a , b , c , d ;
  register double complex sum1 = 0.0 , sum2 = 0.0 ;
  for( abcd = 0 ; abcd < NCNC * NCNC ; abcd++ ) {
    // get the color indices from linearised index abcd
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // [ Uac.Cg5.Dbd.Cg5T.Bac.Cgmu.Bbd.CgmuT ]
    sum1 += c4( U , a , c , Cg5 ,  
		D , b , d , Cg5T ,
		B , a , c , Cgmu , 
		B , b , d , CgmuT ) ;
    // [ Uac.Cg5.Dbd.Cg5T.Bad.Cgmu.Bbc.CgmuT ]
    sum2 += c4( U , a , c , Cg5 ,  
		D , b , d , Cg5T ,
		B , a , d , Cgmu , 
		B , b , c , CgmuT ) ;
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
  struct gamma g5dag = gt_Gdag_gt( g5 , GAMMAS ) ;
  struct gamma gidag = gt_Gdag_gt( gi , GAMMAS ) ;

  // loop all colors
  size_t abcd , a , b , c , d ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    // get the color indices from linearised index abcd
    get_abcd( &a , &b , &c , &d , abcd ) ;
    // [ Bac.g5.Uac.gi.Bbd.g5.Dbd.gi ]
    sum1 += c4( B , a , c , g5 , 
		U , a , c , g5dag ,
		B , b , d , gi , 
		D , b , d , gidag ) ;
    // [ Bad.g5.Uac.g5.Bbc.gi.Dbd.gi ]
    sum2 += c4( B , a , d , g5 , 
		U , a , c , g5dag ,
		B , b , c , gi , 
		D , b , d , gidag ) ;
    // [ Bac.g5.Uad.g5.Bbd.g5.Dbc.gi ]
    sum3 += c4( B , a , c , g5 , 
		U , a , d , g5dag ,
		B , b , d , gi , 
		D , b , c , gidag ) ;
    // [ Bad.g5.Uad.g5.Bbc.gi.Dbc.gi ]
    sum4 += c4( B , a , d , g5 , 
		U , a , d , g5dag ,
		B , b , c , gi , 
		D , b , c , gidag ) ;
  }
  return sum1 - sum2 - sum3 + sum4 ;
}
