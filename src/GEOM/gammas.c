/**
   @file gammas.c
   @brief reshuffler for the gamma matrices

   int gamma_matrix: 
   	gamma matrices for mesons, returns reshuffled index and sign

MAINZ

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
   11-> gamma_0 gamma_2
   12-> gamma_0 gamma_3
   13-> gamma_1 gamma_2
   14-> gamma_1 gamma_3
   15-> gamma_2 gamma_3  

   Mapping of the sign:

   0-> ( 1 + 0 * I)
   1-> ( 0 + 1 * I)
   2-> (-1 + 0 * I)
   3-> ( 0 - 1 * I)

 */
#include "common.h" // needed for struct gamma definition

// gamma multiply
inline void
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
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    const int j = b.ig[i] ; // non-zero column
    a -> ig[i] = c.ig[j] ;
    a -> g[i] = ( b.g[i] + c.g[j] ) & 3 ;
  }
  return ;
#endif
}

// gamma matrix comparison
static int
gamma_comparison( const struct gamma G1 , 
		  const struct gamma G2 ,
		  const char *message )
{
  int j ;
  for( j = 0 ; j < NS ; j++ ) {
    if( G1.ig[j] != G2.ig[j] ||
	G1.g[j] != G2.g[j] ) {
      printf( "%s\n" , message ) ;
      return 1 ;
    }
  }
  return 0 ;
}

// check our gamma matrices are sensible
static int
check_gammas( const struct gamma *GAMMA ) 
{
  // multiply by identity
  struct gamma res ;
  int mu , flag = 0 ;
  for( mu = 0 ; mu < NSNS ; mu++ ) {
    gamma_mmul( &res , GAMMA[ IDENTITY ] , GAMMA[ mu ] ) ;
    flag += gamma_comparison( res , GAMMA[ mu ] , "Identity Multiply" ) ;
    gamma_mmul( &res , GAMMA[ mu ] , GAMMA[ IDENTITY ] ) ;
    flag += gamma_comparison( res , GAMMA[ mu ] , "Identity Multiply" ) ;
  }

  // check gamma_5 * gamma_5 is the identity
  gamma_mmul( &res , GAMMA[ GAMMA_0 ] , GAMMA[ GAMMA_0 ] ) ;
  flag += gamma_comparison( res , GAMMA[ IDENTITY ] , "Gamma0 Gamma0" ) ;
  gamma_mmul( &res , GAMMA[ GAMMA_1 ] , GAMMA[ GAMMA_1 ] ) ;
  flag += gamma_comparison( res , GAMMA[ IDENTITY ] , "Gamma1 Gamma1" ) ;
  gamma_mmul( &res , GAMMA[ GAMMA_2 ] , GAMMA[ GAMMA_2 ] ) ;
  flag += gamma_comparison( res , GAMMA[ IDENTITY ] , "Gamma2 Gamma2" ) ;
  gamma_mmul( &res , GAMMA[ GAMMA_3 ] , GAMMA[ GAMMA_3 ] ) ;
  flag += gamma_comparison( res , GAMMA[ IDENTITY ] , "Gamma3 Gamma3" ) ;
  gamma_mmul( &res , GAMMA[ GAMMA_5 ] , GAMMA[ GAMMA_5 ] ) ;
  flag += gamma_comparison( res , GAMMA[ IDENTITY ] , "Gamma5 Gamma5" ) ;

  // compute gamma_0.gamma_1.gamma_2.gamma_3 = gamma_5 and compare
  struct gamma t1 , t2 ;
  gamma_mmul( &t1 , GAMMA[ GAMMA_0 ] , GAMMA[ GAMMA_1 ] ) ;
  gamma_mmul( &t2 , GAMMA[ GAMMA_2 ] , GAMMA[ GAMMA_3 ] ) ;
  gamma_mmul( &res , t1 , t2 ) ;
  flag += gamma_comparison( res , GAMMA[ GAMMA_5 ] , 
			    "Gamma0.Gamma1.Gamma2.Gamma3 = Gamma5" ) ;

  // if any failure happens we leave in disgust
  if( flag != 0 ) {
    return FAILURE ;
  }
  return SUCCESS ;
}

// conjugate a gamma
static inline uint8_t
gconj( const uint8_t g ) 
{
  return ( ( g & 1 ) ? ( g + 2 ) : g ) & 3 ;
}

// takes a conjugate (and not a dagger)
static const struct gamma
gamma_conj( const struct gamma G )
{
  struct gamma res ;
#if NS == 4 
  res.g[ 0 ] = gconj( G.g[ 0 ] ) ; res.ig[ 0 ] = G.ig[ 0 ] ;
  res.g[ 1 ] = gconj( G.g[ 1 ] ) ; res.ig[ 1 ] = G.ig[ 1 ] ;
  res.g[ 2 ] = gconj( G.g[ 2 ] ) ; res.ig[ 2 ] = G.ig[ 2 ] ; 
  res.g[ 3 ] = gconj( G.g[ 3 ] ) ; res.ig[ 3 ] = G.ig[ 3 ] ;
#else
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    res.ig[ i ] = G.ig[ i ] ;
    res.g[ i ] = gconj( G.g[ i ] ) ;
  }
#endif
  return res ;
}

// computes GAMMA_T . GAMMA_Y . GAMMA_\mu
const struct gamma
CGmu( const struct gamma GAMMA_MU , 
      const struct gamma *GAMMAS )
{
  struct gamma res , tmp ;
  /////////// res  =       gamma_t     .  gamma_y 
  gamma_mmul( &tmp , GAMMAS[ GAMMA_3 ] , GAMMAS[ GAMMA_1 ] ) ;
  gamma_mmul( &res , tmp , GAMMA_MU ) ;
  return res ;
}

// computes ( GAMMA_T . ( C Gamma_\mu )^* . GAMMA_T )
const struct gamma
CGmuT( const struct gamma Cgmu , 
       const struct gamma *GAMMAS )
{
  struct gamma res , tmp ;
  gamma_mmul( &tmp , gamma_conj( Cgmu ) , GAMMAS[ GAMMA_3 ] ) ;
  gamma_mmul( &res , GAMMAS[ GAMMA_3 ] , tmp ) ;
  return res ;
}

#ifdef CHROMA_DIRAC_CONVENTION

// degrand-rossi
void 
make_gammas( struct gamma *GAMMA ,
	     const proptype prop )
{
  // first gamma is the identity
  int s ;
  for( s = 0 ; s < NS ; s++ ) {
    GAMMA[0].ig[s] = s ; // on the diagonal
    GAMMA[0].g[s] = 0 ;  // all ones
  }

  switch( prop ) {
  case NREL :
    // gamma_0 -- x direction
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] =  3 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] =  3 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] =  1 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] =  1 ;
    // gamma_1 -- y direction
    GAMMA[2].ig[0] = 3 ; GAMMA[2].g[0] = 2 ;
    GAMMA[2].ig[1] = 2 ; GAMMA[2].g[1] = 0 ;
    GAMMA[2].ig[2] = 1 ; GAMMA[2].g[2] = 0 ;
    GAMMA[2].ig[3] = 0 ; GAMMA[2].g[3] = 2 ;
    // gamma_2 -- z direction
    GAMMA[4].ig[0] = 2 ; GAMMA[4].g[0] = 3 ;
    GAMMA[4].ig[1] = 3 ; GAMMA[4].g[1] = 1 ;
    GAMMA[4].ig[2] = 0 ; GAMMA[4].g[2] = 1 ;
    GAMMA[4].ig[3] = 1 ; GAMMA[4].g[3] = 3 ;
    // gamma_3 -- t direction
    GAMMA[8].ig[0] = 0 ; GAMMA[8].g[0] = 0 ;
    GAMMA[8].ig[1] = 1 ; GAMMA[8].g[1] = 0 ;
    GAMMA[8].ig[2] = 2 ; GAMMA[8].g[2] = 2 ;
    GAMMA[8].ig[3] = 3 ; GAMMA[8].g[3] = 2 ;
    break ;
  case CHIRAL :
    // gamma_0 is encoded via
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] = 1 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] = 1 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] = 3 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] = 3 ;
    // gamma_1 is encoded via (is this right?)
    GAMMA[2].ig[0] = 3 ; GAMMA[2].g[0] = 2 ;
    GAMMA[2].ig[1] = 2 ; GAMMA[2].g[1] = 0 ;
    GAMMA[2].ig[2] = 1 ; GAMMA[2].g[2] = 0 ;
    GAMMA[2].ig[3] = 0 ; GAMMA[2].g[3] = 2 ;
    // gamma_2 is encoded via
    GAMMA[4].ig[0] = 2 ; GAMMA[4].g[0] = 1 ;
    GAMMA[4].ig[1] = 3 ; GAMMA[4].g[1] = 3 ;
    GAMMA[4].ig[2] = 0 ; GAMMA[4].g[2] = 3 ;
    GAMMA[4].ig[3] = 1 ; GAMMA[4].g[3] = 1 ;
    // gamma_3 is encoded via
    GAMMA[8].ig[0] = 2 ; GAMMA[8].g[0] = 0 ;
    GAMMA[8].ig[1] = 3 ; GAMMA[8].g[1] = 0 ;
    GAMMA[8].ig[2] = 0 ; GAMMA[8].g[2] = 0 ;
    GAMMA[8].ig[3] = 1 ; GAMMA[8].g[3] = 0 ;
    break ;
  }

  // fill in the rest of the chroma basis
  int outer , inner , sub ;
  for( outer = 0 ; outer < NS-1 ; outer++ ) {
    sub = 2 << outer ;
    for( inner = 1 ; inner < sub ; inner++ ) {
      gamma_mmul( &GAMMA[ sub + inner ] , GAMMA[inner] , GAMMA[ sub ] ) ;
    }
  }

  return check_gammas( GAMMA ) ;
}

#else

int
make_gammas( struct gamma *GAMMA ,
	     const proptype prop )
{
  switch( prop ) {
  case NREL :
    // gamma_0 -- x direction
    GAMMA[0].ig[0] = 3 ; GAMMA[0].g[0] =  3 ;
    GAMMA[0].ig[1] = 2 ; GAMMA[0].g[1] =  3 ;
    GAMMA[0].ig[2] = 1 ; GAMMA[0].g[2] =  1 ;
    GAMMA[0].ig[3] = 0 ; GAMMA[0].g[3] =  1 ;
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
  case CHIRAL :
    // gamma_0
    GAMMA[0].ig[0] = 3 ; GAMMA[0].g[0] =  3 ;
    GAMMA[0].ig[1] = 2 ; GAMMA[0].g[1] =  3 ;
    GAMMA[0].ig[2] = 1 ; GAMMA[0].g[2] =  1 ;
    GAMMA[0].ig[3] = 0 ; GAMMA[0].g[3] =  1 ;
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
    // gamma_5 gets flipped no?
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

  // gamma_0 gamma_2 
  gamma_mmul( &GAMMA[11] , GAMMA[0] , GAMMA[2] ) ;

  // gamma_0 gamma_3 
  gamma_mmul( &GAMMA[12] , GAMMA[0] , GAMMA[3] ) ;

  // gamma_1 gamma_2 
  gamma_mmul( &GAMMA[13] , GAMMA[1] , GAMMA[2] ) ;

  // gamma_1 gamma_3 
  gamma_mmul( &GAMMA[14] , GAMMA[1] , GAMMA[3] ) ;

  // gamma_2 gamma_3 
  gamma_mmul( &GAMMA[15] , GAMMA[2] , GAMMA[3] ) ;
  
  return check_gammas( GAMMA ) ;
}

#endif

