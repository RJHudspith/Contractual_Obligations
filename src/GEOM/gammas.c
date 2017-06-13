/**
   @file gammas.c
   @brief gamma matrix business

   0 -> gamma_0 ( x )
   1 -> gamma_1 ( y )
   2 -> gamma_2 ( z )
   3 -> gamma_3 ( time-component )
   4 -> unity
   5 -> gamma_5
   6 -> gamma_0 gamma_5
   7 -> gamma_1 gamma_5
   8 -> gamma_2 gamma_5
   9 -> gamma_3 gamma_5

   10-> gamma_0 gamma_1
   11-> gamma_1 gamma_2
   12-> gamma_2 gamma_0

   13-> gamma_0 gamma_3
   14-> gamma_1 gamma_3
   15-> gamma_2 gamma_3  

   Mapping of the sign:

   0-> ( 1 + 0 * I)
   1-> ( 0 + 1 * I)
   2-> (-1 + 0 * I)
   3-> ( 0 - 1 * I)

 */
#include "common.h" // needed for struct gamma definition

#include "gammas.h"

// computes i GAMMA_Y . GAMMA_T . GAMMA_\mu // gattringer & lang
struct gamma
CGmu( const struct gamma G , 
      const struct gamma *GAMMAS )
{
  struct gamma res , C ;
  // res  =  i .  gamma_y . gamma_t
  gamma_mmul( &C , GAMMAS[ GAMMA_Y ] , GAMMAS[ GAMMA_T ] ) ;
  gamma_muli( &C ) ;
  gamma_mmul( &res , C , G ) ;
  return res ;
}

// takes a conjugate (and not a dagger)
struct gamma
gamma_conj( const struct gamma G )
{
  struct gamma res ;
#if NS == 4 
  res.g[ 0 ] = gconj( G.g[ 0 ] ) ; res.ig[ 0 ] = G.ig[ 0 ] ;
  res.g[ 1 ] = gconj( G.g[ 1 ] ) ; res.ig[ 1 ] = G.ig[ 1 ] ;
  res.g[ 2 ] = gconj( G.g[ 2 ] ) ; res.ig[ 2 ] = G.ig[ 2 ] ; 
  res.g[ 3 ] = gconj( G.g[ 3 ] ) ; res.ig[ 3 ] = G.ig[ 3 ] ;
#else
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    res.ig[ i ] = G.ig[ i ] ;
    res.g[ i ] = gconj( G.g[ i ] ) ;
  }
#endif
  return res ;
}

// takes a conjugate transpose
struct gamma
gamma_dag( const struct gamma G )
{
  struct gamma res ;
#if NS == 4 
  res.g[ 0 ] = gconj( G.g[ G.ig[ 0 ] ] ) ; res.ig[ G.ig[ 0 ] ] = 0 ;
  res.g[ 1 ] = gconj( G.g[ G.ig[ 1 ] ] ) ; res.ig[ G.ig[ 1 ] ] = 1 ;
  res.g[ 2 ] = gconj( G.g[ G.ig[ 2 ] ] ) ; res.ig[ G.ig[ 2 ] ] = 2 ;
  res.g[ 3 ] = gconj( G.g[ G.ig[ 3 ] ] ) ; res.ig[ G.ig[ 3 ] ] = 3 ;
#else
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    res.ig[ G.ig[ i ] ] = i ;
    res.g[ i ] = gconj( G.g[ G.ig[ i ] ] ) ;
  }
#endif
  return res ;
}

// multiply by i
void
gamma_muli( struct gamma *G )
{
  G -> g[0] = ( G -> g[0] + 1 ) & 3 ;
  G -> g[1] = ( G -> g[1] + 1 ) & 3 ;
  G -> g[2] = ( G -> g[2] + 1 ) & 3 ;
  G -> g[3] = ( G -> g[3] + 1 ) & 3 ;
  return ;
}

// multiply by -1
void
gamma_mul_minus1( struct gamma *G )
{
  G -> g[0] = ( G -> g[0] + 2 ) & 3 ;
  G -> g[1] = ( G -> g[1] + 2 ) & 3 ;
  G -> g[2] = ( G -> g[2] + 2 ) & 3 ;
  G -> g[3] = ( G -> g[3] + 2 ) & 3 ;
  return ;
}

// multiply by -i
void
gamma_mul_minusi( struct gamma *G )
{
  G -> g[0] = ( G -> g[0] + 3 ) & 3 ;
  G -> g[1] = ( G -> g[1] + 3 ) & 3 ;
  G -> g[2] = ( G -> g[2] + 3 ) & 3 ;
  G -> g[3] = ( G -> g[3] + 3 ) & 3 ;
  return ;
}

// gamma multiply
void
gamma_mmul( struct gamma *__restrict a ,
	    const struct gamma b ,
	    const struct gamma c )
{
#if NS == 4
  a -> ig[0] = c.ig[ b.ig[0] ] ;
  a ->  g[0] = ( b.g[0] + c.g[ b.ig[0] ] ) & 3 ;
  a -> ig[1] = c.ig[ b.ig[1] ] ;
  a ->  g[1] = ( b.g[1] + c.g[ b.ig[1] ] ) & 3 ;
  a -> ig[2] = c.ig[ b.ig[2] ] ;
  a ->  g[2] = ( b.g[2] + c.g[ b.ig[2] ] ) & 3 ;
  a -> ig[3] = c.ig[ b.ig[3] ] ;
  a ->  g[3] = ( b.g[3] + c.g[ b.ig[3] ] ) & 3 ;
#else
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    const size_t j = b.ig[i] ; // non-zero column
    a -> ig[i] = c.ig[j] ;
    a -> g[i] = ( b.g[i] + c.g[j] ) & 3 ;
  }
  return ;
#endif
}

// takes a transpose
struct gamma
gamma_transpose( const struct gamma G )
{
  struct gamma res ;
#if NS == 4 
  res.g[ 0 ] = G.g[ G.ig[ 0 ] ] ; res.ig[ G.ig[ 0 ] ] = 0 ;
  res.g[ 1 ] = G.g[ G.ig[ 1 ] ] ; res.ig[ G.ig[ 1 ] ] = 1 ;
  res.g[ 2 ] = G.g[ G.ig[ 2 ] ] ; res.ig[ G.ig[ 2 ] ] = 2 ;
  res.g[ 3 ] = G.g[ G.ig[ 3 ] ] ; res.ig[ G.ig[ 3 ] ] = 3 ;
#else
  size_t i ;
  for( i = 0 ; i < NS ; i++ ) {
    res.ig[ G.ig[ i ] ] = i ;
    res.g[ i ] = G.g[ G.ig[ i ] ] ;
  }
#endif
  return res ;
}

// conjugate a gamma
uint8_t
gconj( const uint8_t g ) 
{
  return ( ( g % 2 ) == 1 ? ( g + 2 ) : g ) & 3 ;
}

// computes ( GAMMA_T . ( G )^* . GAMMA_T )
struct gamma
gt_Gconj_gt( const struct gamma G, 
	     const struct gamma gt )
{
  struct gamma res , tmp ;
  gamma_mmul( &tmp , gamma_conj( G ) , gt ) ;
  gamma_mmul( &res , gt , tmp ) ;
  return res ;
}

// computes ( GAMMA_T . ( G )^\dagger . GAMMA_T )
struct gamma
gt_Gdag_gt( const struct gamma G , 
	    const struct gamma gt )
{
  struct gamma res , tmp ;
  gamma_mmul( &tmp , gamma_dag( G ) , gt) ;
  gamma_mmul( &res , gt, tmp ) ;
  return res ;
}

// create a look up table for our gamma matrices
int
make_gammas( struct gamma *GAMMA ,
	     const proptype prop )
{
  switch( prop ) {
  case NREL_FWD :
  case NREL_BWD :
    // gamma_0 -- x direction
    GAMMA[0].ig[0] = 3 ; GAMMA[0].g[0] = 3 ;
    GAMMA[0].ig[1] = 2 ; GAMMA[0].g[1] = 3 ;
    GAMMA[0].ig[2] = 1 ; GAMMA[0].g[2] = 1 ;
    GAMMA[0].ig[3] = 0 ; GAMMA[0].g[3] = 1 ;
    // gamma_1 -- y direction
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] = 2 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] = 0 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] = 0 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] = 2 ;
    // gamma_2 -- z direction
    GAMMA[2].ig[0] = 2 ; GAMMA[2].g[0] = 3 ;
    GAMMA[2].ig[1] = 3 ; GAMMA[2].g[1] = 1 ;
    GAMMA[2].ig[2] = 0 ; GAMMA[2].g[2] = 1 ;
    GAMMA[2].ig[3] = 1 ; GAMMA[2].g[3] = 3 ;
    // gamma_3 -- t direction
    GAMMA[3].ig[0] = 0 ; GAMMA[3].g[0] = 0 ;
    GAMMA[3].ig[1] = 1 ; GAMMA[3].g[1] = 0 ;
    GAMMA[3].ig[2] = 2 ; GAMMA[3].g[2] = 2 ;
    GAMMA[3].ig[3] = 3 ; GAMMA[3].g[3] = 2 ;
    // gamma_5 
    GAMMA[5].ig[0] = 2 ; GAMMA[5].g[0] = 2 ;
    GAMMA[5].ig[1] = 3 ; GAMMA[5].g[1] = 2 ;
    GAMMA[5].ig[2] = 0 ; GAMMA[5].g[2] = 2 ;
    GAMMA[5].ig[3] = 1 ; GAMMA[5].g[3] = 2 ;
    break ;
  case STATIC :
  case CHIRAL :
    // gamma_0
    GAMMA[0].ig[0] = 3 ; GAMMA[0].g[0] = 3 ;
    GAMMA[0].ig[1] = 2 ; GAMMA[0].g[1] = 3 ;
    GAMMA[0].ig[2] = 1 ; GAMMA[0].g[2] = 1 ;
    GAMMA[0].ig[3] = 0 ; GAMMA[0].g[3] = 1 ;
    // gamma_1
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] = 2 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] = 0 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] = 0 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] = 2 ;
    // gamma_2
    GAMMA[2].ig[0] = 2 ; GAMMA[2].g[0] = 3 ;
    GAMMA[2].ig[1] = 3 ; GAMMA[2].g[1] = 1 ;
    GAMMA[2].ig[2] = 0 ; GAMMA[2].g[2] = 1 ;
    GAMMA[2].ig[3] = 1 ; GAMMA[2].g[3] = 3 ;
    // gamma_3
    GAMMA[3].ig[0] = 2 ; GAMMA[3].g[0] = 2 ;
    GAMMA[3].ig[1] = 3 ; GAMMA[3].g[1] = 2 ;
    GAMMA[3].ig[2] = 0 ; GAMMA[3].g[2] = 2 ;
    GAMMA[3].ig[3] = 1 ; GAMMA[3].g[3] = 2 ;
    // gamma_5
    GAMMA[5].ig[0] = 0 ; GAMMA[5].g[0] = 2 ;
    GAMMA[5].ig[1] = 1 ; GAMMA[5].g[1] = 2 ;
    GAMMA[5].ig[2] = 2 ; GAMMA[5].g[2] = 0 ;
    GAMMA[5].ig[3] = 3 ; GAMMA[5].g[3] = 0 ;
    break ;
  }

  // unity 
  GAMMA[4].ig[0] = 0 ; GAMMA[4].g[0] = 0 ;
  GAMMA[4].ig[1] = 1 ; GAMMA[4].g[1] = 0 ;
  GAMMA[4].ig[2] = 2 ; GAMMA[4].g[2] = 0 ;
  GAMMA[4].ig[3] = 3 ; GAMMA[4].g[3] = 0 ;

  // and multiply out the rest
  // gamma_0 gamma_5 
  gamma_mmul( &GAMMA[6] , GAMMA[0] , GAMMA[5] ) ;

  // gamma_1 gamma_5 
  gamma_mmul( &GAMMA[7] , GAMMA[1] , GAMMA[5] ) ;

  // gamma_2 gamma_5 
  gamma_mmul( &GAMMA[8] , GAMMA[2] , GAMMA[5] ) ;

  // gamma_3 gamma_5 
  gamma_mmul( &GAMMA[9] , GAMMA[3] , GAMMA[5] ) ;

  // gamma_0 gamma_1 
  gamma_mmul( &GAMMA[10] , GAMMA[0] , GAMMA[1] ) ;

  // gamma_1 gamma_2 
  gamma_mmul( &GAMMA[11] , GAMMA[1] , GAMMA[2] ) ;

  // gamma_2 gamma_0
  gamma_mmul( &GAMMA[12] , GAMMA[2] , GAMMA[0] ) ;

  // gamma_0 gamma_3 
  gamma_mmul( &GAMMA[13] , GAMMA[0] , GAMMA[3] ) ;

  // gamma_1 gamma_3 
  gamma_mmul( &GAMMA[14] , GAMMA[1] , GAMMA[3] ) ;

  // gamma_2 gamma_3 
  gamma_mmul( &GAMMA[15] , GAMMA[2] , GAMMA[3] ) ;
  
  // compliance is now checked in gamma_tests.c as part 
  // of the unit-testing structure
  return SUCCESS ;
}

// have a look at our gammas
void
picture_gamma( const struct gamma G ) 
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      // 
      if( G.ig[ d1 ] == d2 ) {
	switch( G.g[d1] ) {
	case 0 : fprintf( stdout , " +1 " ) ; break ;
	case 1 : fprintf( stdout , " +i " ) ; break ;
	case 2 : fprintf( stdout , " -1 " ) ; break ;
	case 3 : fprintf( stdout , " -i " ) ; break ;
	}
      } else {
	fprintf( stdout , "  0 " ) ;
      }
      //
    }
    fprintf( stdout , "\n" ) ;
  }
  fprintf( stdout , "\n" ) ;
  return ;
}

// setup our gamma matrices
int
setup_gamma( struct gamma *GAMMAS , 
	     const struct propagator *prop , 
	     const size_t Nprops )
{
  // catch this nonsense
  if( Nprops == 0 ) {
    return FAILURE ;
  }

  // fall through if chiral
  proptype basis = CHIRAL ;

  size_t mu ;
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    // if we hit a single NREL prop we set gamma basis to NREL
    if( prop[mu].basis == NREL_FWD || prop[mu].basis == NREL_BWD ) {
      basis = NREL_FWD ;
      break ;
    }
  }

  return make_gammas( GAMMAS , basis ) ;
}
