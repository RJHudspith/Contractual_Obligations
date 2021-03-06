/**
   @file dibaryon_contractions.c
   @brief SU(2) dibaryon contraction code
 */
#include "common.h"

#include "gammas.h"             
#include "Ospinor.h"
#include "spinmatrix_ops.h" 
#include "spinor_ops.h"     // transpose_spinor()

// does Tr[ S1T G1 S2 G2 S3T G3 S4 G4 ]
static inline double complex
single_tr( struct spinmatrix block1 ,
	   struct spinmatrix block2 )
{
  return trace_prod_spinmatrices( (const void*)block1.D ,
				  (const void*)block2.D ) ;
}

// does Tr[ S1T G1 S2 G2 ] x Tr[ S3T G3 S4 G4 ]
static inline double complex
double_tr( struct spinmatrix block1 ,
	   struct spinmatrix block2 )
{
  return
    spinmatrix_trace( (const void*)block1.D ) *
    spinmatrix_trace( (const void*)block2.D ) ;
}

// from lexicographical index abcd get index a,b,c,d --> only correct for NC=2
static inline void
get_abcd( size_t *a , size_t *b , size_t *ap , size_t *bp , const size_t abcd )
{
  *a  = abcd >> 3 ;
  *b  = ( abcd >> 2 )&1 ;
  *ap = ( abcd >> 1 )&1 ;
  *bp = abcd&1 ;
  return ;
}

// given a,b,ap,bp return lexicographical index abcd
static inline size_t
get_idx( const size_t a , const size_t b , const size_t ap , const size_t bp )
{
  return bp + NC*( ap + NC*( b + NC*a ) ) ;
}

double complex
dibaryon_contract( struct spinor S1 ,
		   struct spinor S2 ,
		   const struct gamma *GAMMAS ,
		   const size_t GSRC ,
		   const size_t GSNK )
{
  // precomputations -> swap color and dirac indices to expose spinmatrices
  struct Ospinor OS1T = spinor_to_Ospinor( transpose_spinor( S1 ) ) ;
  struct Ospinor OS1  = spinor_to_Ospinor( S1 ) ;
  struct Ospinor OS2T = spinor_to_Ospinor( transpose_spinor( S2 ) ) ;
  struct Ospinor OS2  = spinor_to_Ospinor( S2 ) ;
  
  // gamma precomputations
  const struct gamma G1   = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
  const struct gamma tG1t = gt_Gdag_gt( G1 , GAMMAS[ GAMMA_T ] ) ;

  // gamma precomputations
  const struct gamma G2   = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
  const struct gamma tG2t = gt_Gdag_gt( G2 , GAMMAS[ GAMMA_T ] ) ;
  
  // precompute the diquark blocks, called blk
  struct spinmatrix blk11[ NCNC*NCNC ] , blk12[ NCNC*NCNC ] ,
    blk21[ NCNC*NCNC ] , blk22[ NCNC*NCNC ] ;

  // multiply left and right by gammas into tmp with color indices ab
  // and multiply on the left by the transposed prop and save result
  struct spinmatrix tmp ;
  size_t a , b , ap , bp , abcd ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &ap , &bp , abcd ) ;
    // g1g1 case
    tmp = OS2.C[a][b] ;
    gamma_spinmatrix_lr( &tmp , G1 , tG1t ) ;
    spinmatrix_multiply( (void*)blk11[abcd].D ,
			 (void*)OS2T.C[ap][bp].D , (void*)tmp.D ) ;
    // g1g2 case
    tmp = OS2.C[a][b] ;
    gamma_spinmatrix_lr( &tmp , G1 , tG2t ) ;
    spinmatrix_multiply( (void*)blk12[abcd].D ,
			 (void*)OS2T.C[ap][bp].D , (void*)tmp.D ) ;
    // g2g1 case
    tmp = OS1.C[a][b] ;
    gamma_spinmatrix_lr( &tmp , G2 , tG1t ) ;
    spinmatrix_multiply( (void*)blk21[abcd].D ,
			 (void*)OS1T.C[ap][bp].D , (void*)tmp.D ) ;
    // g2g2 case
    tmp = OS1.C[a][b] ;
    gamma_spinmatrix_lr( &tmp , G2 , tG2t ) ;
    spinmatrix_multiply( (void*)blk22[abcd].D ,
			 (void*)OS1T.C[ap][bp].D , (void*)tmp.D ) ;
  }

  // do all of the 4! contractions summing over the 16 color combinations
  // note that I use the identity that G^T for gamma_i is G
  register double complex sum = 0.0 ;
  for( abcd = 0 ; abcd < NCNC*NCNC ; abcd++ ) {
    get_abcd( &a , &b , &ap , &bp , abcd ) ;
    // diagram 1/24
    sum += double_tr( blk11[get_idx(a,ap,b,bp)] , blk22[get_idx(a,ap,b,bp)] ) ;
    // diagram 2/24
    sum -= double_tr( blk11[get_idx(b,ap,a,bp)] , blk22[get_idx(a,ap,b,bp)] ) ;
    // diagram 3/24
    sum += single_tr( blk11[get_idx(b,ap,a,ap)] , blk22[get_idx(a,bp,b,bp)] ) ;
    // diagram 4/24
    sum -= single_tr( blk11[get_idx(a,ap,b,ap)] , blk22[get_idx(a,bp,b,bp)] ) ;
    // diagram 5/24
    sum -= single_tr( blk11[get_idx(b,ap,a,bp)] , blk22[get_idx(a,ap,b,bp)] ) ;
    // diagram 6/24
    sum += single_tr( blk11[get_idx(a,ap,b,bp)] , blk22[get_idx(a,ap,b,bp)] ) ;
    // diagram 7/24
    sum -= double_tr( blk11[get_idx(a,ap,b,bp)] , blk22[get_idx(b,ap,a,bp)] ) ;
    // diagram 8/24
    sum += double_tr( blk11[get_idx(b,ap,a,bp)] , blk22[get_idx(b,ap,a,bp)] ) ;
    // diagram 9/24
    sum -= single_tr( blk12[get_idx(a,ap,b,bp)] , blk21[get_idx(b,ap,a,bp)] ) ;
    // diagram 10/24
    sum += single_tr( blk12[get_idx(b,ap,a,bp)] , blk21[get_idx(b,ap,a,bp)] ) ;
    // diagram 11/24
    sum += single_tr( blk11[get_idx(b,bp,a,bp)] , blk22[get_idx(a,ap,b,ap)] ) ;
    // diagram 12/24
    sum -= single_tr( blk11[get_idx(a,bp,b,bp)] , blk22[get_idx(a,ap,b,ap)] ) ;
    //diagram 13/24
    sum -= single_tr( blk12[get_idx(a,ap,b,ap)] , blk21[get_idx(a,bp,b,bp)] ) ;
    // diagram 14/24
    sum += single_tr( blk11[get_idx(a,ap,b,ap)] , blk22[get_idx(b,bp,a,bp)] ) ;
    // diagram 15/24
    sum += single_tr( blk12[get_idx(a,ap,b,bp)] , blk21[get_idx(a,ap,b,bp)] ) ;
    // diagram 16/24
    sum -= single_tr( blk12[get_idx(b,ap,a,bp)] , blk21[get_idx(a,ap,b,bp)] ) ;
    // diagram 17/24
    sum += double_tr( blk12[get_idx(a,ap,b,bp)] , blk21[get_idx(a,ap,b,bp)] ) ;
    // diagram 18/24
    sum -= double_tr( blk12[get_idx(b,ap,a,bp)] , blk21[get_idx(a,ap,b,bp)] ) ;
    // diagram 19/24
    sum += single_tr( blk11[get_idx(b,ap,a,bp)] , blk22[get_idx(b,ap,a,bp)] ) ;
    // diagram 20/24
    sum -= single_tr( blk11[get_idx(a,ap,b,bp)] , blk22[get_idx(b,ap,a,bp)] ) ;
    // diagram 21/24
    sum -= single_tr( blk12[get_idx(a,bp,b,bp)] , blk21[get_idx(a,ap,b,ap)] ) ;
    // diagram 22/24
    sum += single_tr( blk12[get_idx(b,bp,a,bp)] , blk21[get_idx(a,ap,b,ap)] ) ;
    // diagram 23/24
    sum -= double_tr( blk12[get_idx(a,ap,b,bp)] , blk21[get_idx(b,ap,a,bp)] ) ;
    // diagram 24/24
    sum += double_tr( blk12[get_idx(b,ap,a,bp)] , blk21[get_idx(b,ap,a,bp)] ) ;
  }
  
  return sum ;
}

// slower version down here with explicit gamma stuff -- used for checking
#if 0

// does Tr[ S1T G1 S2 G2 S3T G3 S4 G4 ]
double complex
single_tr2( const struct spinmatrix S1T ,
	   const struct gamma G1 ,
	   const struct spinmatrix S2 ,
	   const struct gamma G2 ,
	   const struct spinmatrix S3T ,
	   const struct gamma G3 ,
	   const struct spinmatrix S4 ,
	   const struct gamma G4 )
{
  struct spinmatrix B = S2 , D = S4 , E , F ;
  gamma_spinmatrix_lr( &B , G1 , G2 ) ;
  spinmatrix_multiply( (void*)E.D , (void*)S1T.D , (void*)B.D ) ;
  gamma_spinmatrix_lr( &D , G3 , G4 ) ;
  spinmatrix_multiply( (void*)F.D , (void*)S3T.D , (void*)D.D ) ;
  return trace_prod_spinmatrices( (const void*)E.D , (const void*)F.D ) ;
}

// does Tr[ S1T G1 S2 G2 ] x Tr[ S3T G3 S4 G4 ]
double complex
double_tr2( const struct spinmatrix S1T ,
	   const struct gamma G1 ,
	   const struct spinmatrix S2 ,
	   const struct gamma G2 ,
	   const struct spinmatrix S3T ,
	   const struct gamma G3 ,
	   const struct spinmatrix S4 ,
	   const struct gamma G4 )
{
  struct spinmatrix A = S1T , B = S2 , C = S3T , D = S4 ;
  gamma_spinmatrix_lr( &B , G1 , G2 ) ;
  gamma_spinmatrix_lr( &D , G3 , G4 ) ;
  return  trace_prod_spinmatrices( (const void*)A.D , (const void*)B.D )
    * trace_prod_spinmatrices( (const void*)C.D , (const void*)D.D )
    ;
}


double complex
dibaryon_contract2( struct spinor S ,
		    const struct gamma *GAMMAS ,
		    const size_t GSRC ,
		    const size_t GSNK )
{
  // precomputations -> swap color and dirac indices to expose spinmatrices
  struct Ospinor OST = spinor_to_Ospinor( transpose_spinor( S ) ) ;
  struct Ospinor OS  = spinor_to_Ospinor( S ) ;

  // gamma precomputations
  const struct gamma G1    = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
  const struct gamma G1T   = gamma_transpose( G1 ) ;
  const struct gamma tG1t  = gt_Gdag_gt( G1 , GAMMAS[ GAMMA_T ] ) ;
  const struct gamma tG1tT = gamma_transpose( tG1t ) ;

  // gamma precomputations
  const struct gamma G2    = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
  const struct gamma G2T   = gamma_transpose( G2 ) ;
  const struct gamma tG2t  = gt_Gdag_gt( G2 , GAMMAS[ GAMMA_T ] ) ;
  const struct gamma tG2tT = gamma_transpose( tG2t ) ;
  
  size_t a , b , ap , bp ;

  double complex sum = 0.0 ;
  for( a = 0 ; a < NC ; a++ ) {
    for( b = 0 ; b < NC ; b++ ) {
      for( ap = 0 ; ap < NC ; ap++ ) {
	for( bp = 0 ; bp < NC ; bp++ ) {

	  // diagram 1/24
	  sum += double_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 2/24
	  sum -= double_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 3/24
	  sum += single_tr2( OST.C[b][ap] , G1T , OS.C[a][ap] , tG1tT ,
			     OST.C[a][bp] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 4/24
	  sum -= single_tr2( OST.C[a][ap] , G1 , OS.C[b][ap] , tG1tT ,
			     OST.C[a][bp] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 5/24
	  sum -= single_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 6/24
	  sum += single_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG2t ) ;
	  // diagram 7/24
	  sum -= double_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG1t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG2t ) ;
	  // diagram 8/24
	  sum += double_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG1t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG2t ) ;
	  // diagram 9/24
	  sum -= single_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG2t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG1t ) ;
	  // diagram 10/24
	  sum += single_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG2t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG1t ) ;
	  // diagram 11/24
	  sum += single_tr2( OST.C[b][bp] , G1T , OS.C[a][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][ap] , tG2tT ) ;
	  // diagram 12/24
	  sum -= single_tr2( OST.C[a][bp] , G1 , OS.C[b][bp] , tG1t ,
			     OST.C[a][ap] , G2 , OS.C[b][ap] , tG2tT ) ;
	  //diagram 13/24
	  sum -= single_tr2( OST.C[a][ap] , G1 , OS.C[b][ap] , tG2tT ,
			     OST.C[a][bp] , G2 , OS.C[b][bp] , tG1t ) ;
	  // diagram 14/24
	  sum += single_tr2( OST.C[a][ap] , G1 , OS.C[b][ap] , tG1tT ,
			     OST.C[b][bp] , G2T , OS.C[a][bp] , tG2t ) ;
	  // diagram 15/24
	  sum += single_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG1t ) ;
	  // diagram 16/24
	  sum -= single_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG1t ) ;
	  // diagram 17/24
	  sum += double_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG1t ) ;
	  // diagram 18/24
	  sum -= double_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][bp] , tG1t ) ;
	  // diagram 19/24
	  sum += single_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG1t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG2t ) ;
	  // diagram 20/24
	  sum -= single_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG1t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG2t ) ;
	  // diagram 21/24
	  sum -= single_tr2( OST.C[a][bp] , G1 , OS.C[b][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][ap] , tG1tT ) ;
	  // diagram 22/24
	  sum += single_tr2( OST.C[b][bp] , G1T , OS.C[a][bp] , tG2t ,
			     OST.C[a][ap] , G2 , OS.C[b][ap] , tG1tT ) ;
	  // diagram 23/24
 	  sum -= double_tr2( OST.C[a][ap] , G1 , OS.C[b][bp] , tG2t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG1t ) ;
	  // diagram 24/24
	  sum += double_tr2( OST.C[b][ap] , G1T , OS.C[a][bp] , tG2t ,
			     OST.C[b][ap] , G2T , OS.C[a][bp] , tG1t ) ;
	}
      }
    }
  }
  
  return sum ;
}

#endif
