/**
   @file brutal_mesons.c
   @brief brute-force meson contractions
 */
#include "common.h"

#include "mesons.h"      // for the allocate and free correlator functions
#include "spinor_ops.h"  // bilinear trace

double complex
meson_trace2( const struct gamma GSNK ,
	      const struct spinor adj ,
	      const struct gamma GSRC ,
	      const struct spinor fwd )
{
  // loop counters
  int d1 , d2 , c1 , c2 ;
  int sign ; // permutation of the fourth roots of unity

  register double complex corr = 0.0 ;
  register double rloc_corr , iloc_corr ;

  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {

      // map the gamma multiplied indices for legibility 
      const int id1 = GSRC.ig[ d1 ] ;
      const int id2 = GSNK.ig[ d2 ] ;

      // adjust the sign is a permutation of the fourth roots of unity for all sensible
      // gamma matrix conventions
      sign = ( GSRC.g[d1] + GSNK.g[d2] ) & 3 ;

      rloc_corr = 0.0 ;
      iloc_corr = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  // this does color trace of [ fwd[c1][c2] * adj[c2][c1] ]
	  // real sum
	  rloc_corr += ( creal( fwd.D[d1][d2].C[c1][c2] ) * creal( adj.D[id1][id2].C[c1][c2] ) + \
			 cimag( fwd.D[d1][d2].C[c1][c2] ) * cimag( adj.D[id1][id2].C[c1][c2] ) ) ;
	  //
	  iloc_corr += ( creal( fwd.D[d1][d2].C[c1][c2] ) * cimag( adj.D[id1][id2].C[c1][c2] ) - \
			 cimag( fwd.D[d1][d2].C[c1][c2] ) * creal( adj.D[id1][id2].C[c1][c2] ) ) ;
	}
      }
      // is just a permutation. With this basis this requires NS = 4
      switch( sign ) {
      case 0 : corr += rloc_corr + I * iloc_corr ; break ;
      case 1 : corr += -iloc_corr + I * rloc_corr ; break ;
      case 2 : corr += -rloc_corr - I * iloc_corr ; break ;
      case 3 : corr += iloc_corr - I * rloc_corr ; break ;
      }
    } 
  }
  return -corr ;
}

// dirty computation of mesons
double complex
meson_trace( const struct gamma GSNK ,
	     const struct spinor adj ,
	     const struct gamma GSRC ,
	     const struct spinor fwd )
{
  struct spinor tmp1 = adj ; 
  gamma_mul_l( &tmp1 , GSNK ) ;

  struct spinor tmp2 = fwd ; 
  gamma_mul_l( &tmp2 , GSRC ) ;

  return bilinear_trace( tmp1 , tmp2 ) ;
}

// computes meson correlators
int
single_mesons_bruteforce( FILE *prop1 , 
			  const proptype proptype1 )
{
  // data structure for holding the contractions
  struct correlator **corr = calloc( NS*NS , sizeof( struct correlator* ) ) ;

  allocate_corrs( corr ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1adj = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NS * NS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS , proptype1 ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop1 , S1 , proptype1 ) ;

    // loop spatial hypercube
    int site ;
    #pragma omp parallel for private( site )
    for( site = 0 ; site < VOL3 ; site++ ) {
      // pre-compute the adjoint
      full_adj( &S1adj[site] , S1[ site ] , GAMMAS[ GAMMA_5 ] ) ;
    }

    int GSRC = 0 ;
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NS*NS ; GSRC++ ) {
      
      int GSNK ;
      for( GSNK = 0 ; GSNK < NS*NS ; GSNK++ ) {
	
	register double complex sum = 0.0 ;
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	   sum += meson_trace2( GAMMAS[ GSNK ] , S1adj[site] , GAMMAS[ GSRC ] , S1[ site ] ) ;		      
	}
	corr[ GSRC ][ GSNK ].C[ t ] = (double complex)sum ;
      }
    }
    
    printf("\rContractions done %.f %%",(t+1)/((L0)/100.));fflush(stdout);
  }
  printf("\n");	

  // & do something with the computed correlators
#ifdef DEBUG
  debug_mesons( "LL-mesons" , (const struct correlator**)corr ) ;
#endif

#if 0
  int GS , GN ;
  for( GS = 0 ; GS < NS*NS ; GS++ ){
    for( GN = 0 ; GN < NS*NS ; GN++ ) {

      int t ;
      for( t = 0 ; t < L0 ; t++ ) {
	printf( "[ %d %d %d ] %e %e \n" , GS , GN , t , creal( corr[ GS ][ GN ].C[t] ) , cimag( corr[GS][GN].C[t] )) ;
      }
      //
      printf( "\n" ) ;
    }
  }
#endif

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;
  free( S1adj ) ;

  return SUCCESS ;
}
