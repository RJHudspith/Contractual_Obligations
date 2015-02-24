/**
   @file gammas.c
   @brief reshuffler for the gamma matrices

   int gamma_matrix: 
   	gamma matrices for mesons, returns reshuffled index and sign
 
       0 -> gamma_0 (time-component)
       1 -> gamma_1
       2 -> gamma_2
       3 -> gamma_3
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
static void
gamma_mmul( struct gamma *a ,
	    const struct gamma b ,
	    const struct gamma c )
{
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    const int j = b.ig[i] ; // non-zero column
    a -> ig[i] = c.ig[j] ;
    a -> g[i] = ( b.g[i] + c.g[j] ) & 3 ;
  }
  return ;
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

  // gamma_0 is encoded via
  GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] = 1 ;
  GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] = 1 ;
  GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] = 3 ;
  GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] = 3 ;

  // gamma_1 is encoded via
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
  
  // fill in the rest of the chroma basis
  int outer , inner , sub ;
  for( outer = 0 ; outer < NS-1 ; outer++ ) {
    sub = 2 << outer ;
    for( inner = 1 ; inner < sub ; inner++ ) {
      gamma_mmul( &GAMMA[ sub + inner ] , GAMMA[inner] , GAMMA[ sub ] ) ;
    }
  }
  return ;
}

#else

void 
make_gammas( struct gamma *GAMMA ,
	     const proptype prop )
{
  switch( prop ) {
  case NREL :
  case CHIRAL_TO_NREL :
    // gamma_0
    GAMMA[0].ig[0] = 0 ; GAMMA[0].g[0] = 0 ;
    GAMMA[0].ig[1] = 1 ; GAMMA[0].g[1] = 0 ;
    GAMMA[0].ig[2] = 2 ; GAMMA[0].g[2] = 2 ;
    GAMMA[0].ig[3] = 3 ; GAMMA[0].g[3] = 2 ;

    // gamma_1
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] =  3 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] =  3 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] =  1 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] =  1 ;

    // gamma_2
    GAMMA[2].ig[0] = 3 ; GAMMA[2].g[0] = 2 ;
    GAMMA[2].ig[1] = 2 ; GAMMA[2].g[1] = 0 ;
    GAMMA[2].ig[2] = 1 ; GAMMA[2].g[2] = 0 ;
    GAMMA[2].ig[3] = 0 ; GAMMA[2].g[3] = 2 ;

    // gamma_3
    GAMMA[3].ig[0] = 2 ; GAMMA[3].g[0] = 3 ;
    GAMMA[3].ig[1] = 3 ; GAMMA[3].g[1] = 1 ;
    GAMMA[3].ig[2] = 0 ; GAMMA[3].g[2] = 1 ;
    GAMMA[3].ig[3] = 1 ; GAMMA[3].g[3] = 3 ;
  
    // gamma_5 
    GAMMA[5].ig[0] = 2 ; GAMMA[5].g[0] = 0 ;
    GAMMA[5].ig[1] = 3 ; GAMMA[5].g[1] = 0 ;
    GAMMA[5].ig[2] = 0 ; GAMMA[5].g[2] = 0 ;
    GAMMA[5].ig[3] = 1 ; GAMMA[5].g[3] = 0 ;
    break ;
  case CHIRAL :
    // gamma_0
    GAMMA[0].ig[0] = 2 ; GAMMA[0].g[0] = 2 ;
    GAMMA[0].ig[1] = 3 ; GAMMA[0].g[1] = 2 ;
    GAMMA[0].ig[2] = 0 ; GAMMA[0].g[2] = 2 ;
    GAMMA[0].ig[3] = 1 ; GAMMA[0].g[3] = 2 ;

    // gamma_1
    GAMMA[1].ig[0] = 3 ; GAMMA[1].g[0] =  3 ;
    GAMMA[1].ig[1] = 2 ; GAMMA[1].g[1] =  3 ;
    GAMMA[1].ig[2] = 1 ; GAMMA[1].g[2] =  1 ;
    GAMMA[1].ig[3] = 0 ; GAMMA[1].g[3] =  1 ;

    // gamma_2
    GAMMA[2].ig[0] = 3 ; GAMMA[2].g[0] = 2 ;
    GAMMA[2].ig[1] = 2 ; GAMMA[2].g[1] = 0 ;
    GAMMA[2].ig[2] = 1 ; GAMMA[2].g[2] = 0 ;
    GAMMA[2].ig[3] = 0 ; GAMMA[2].g[3] = 2 ;

    // gamma_3
    GAMMA[3].ig[0] = 2 ; GAMMA[3].g[0] = 3 ;
    GAMMA[3].ig[1] = 3 ; GAMMA[3].g[1] = 1 ;
    GAMMA[3].ig[2] = 0 ; GAMMA[3].g[2] = 1 ;
    GAMMA[3].ig[3] = 1 ; GAMMA[3].g[3] = 3 ;
  
    // gamma_5 
    GAMMA[5].ig[0] = 0 ; GAMMA[5].g[0] = 0 ;
    GAMMA[5].ig[1] = 1 ; GAMMA[5].g[1] = 0 ;
    GAMMA[5].ig[2] = 2 ; GAMMA[5].g[2] = 2 ;
    GAMMA[5].ig[3] = 3 ; GAMMA[5].g[3] = 2 ;
    break ;
  }

  // unity 
  GAMMA[4].ig[0] = 0 ; GAMMA[4].g[0] = 0 ;
  GAMMA[4].ig[1] = 1 ; GAMMA[4].g[1] = 0 ;
  GAMMA[4].ig[2] = 2 ; GAMMA[4].g[2] = 0 ;
  GAMMA[4].ig[3] = 3 ; GAMMA[4].g[3] = 0 ;

  // and multiply out the rest
  // amma_0 gamma_5 
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
  
  return ;
}

#endif

