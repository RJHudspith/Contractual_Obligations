/**
   @file brutal_mesons.c
   @brief brute-force meson contractions
 
   TODO :: make this equivalent to the mesons code
*/
#include "common.h"

#include "contractions.h" // bilinear trace and gamma multiplies
#include "correlators.h"  // for the allocate and free correlator functions
#include "gammas.h"       // make gammas
#include "GLU_timer.h"    // print_time() function
#include "io.h"           // read_prop

// meson trace
double complex
meson_trace( const struct gamma GSNK ,
	     const struct spinor S2 ,
	     const struct gamma GSRC ,
	     const struct spinor S1 )
{
  /*
  struct spinor tmp1 = S2 ;
  gamma_mul_l( &tmp1 , GSNK ) ;

  struct spinor tmp2 = S1 ;
  gamma_mul_l( &tmp2 , GSRC ) ;
  
  return bilinear_trace( tmp1 , tmp2 ) ;
  */
  return simple_meson_contract( GSNK , S2 , GSRC , S1 ) ;
}

// computes meson correlators
int
single_mesons_bruteforce( struct propagator prop ,
			  const char *outfile )
{
  // data structure for holding the contractions
  struct correlator **corr = allocate_corrs( NSNS , NSNS ) ;

  // and our spinor
  struct spinor *S1 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *S1adj = calloc( VOL3 , sizeof( struct spinor ) ) ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  make_gammas( GAMMAS , prop.basis ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    read_prop( prop , S1 ) ;

    // compute the full adjoint
    int i ;
#pragma omp parallel for private(i)
    for( i = 0 ; i < VOL3 ; i++ ) {
      full_adj( &S1adj[ i ] , S1[ i ] , GAMMAS[ GAMMA_5 ] ) ;
    }

    // loop over all source and sink indices
    int GSRC = 0 ;
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NS*NS ; GSRC++ ) {
      
      int GSNK ;
      for( GSNK = 0 ; GSNK < NS*NS ; GSNK++ ) {
	
	register double complex sum = 0.0 ;
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_trace( GAMMAS[ GSNK ] , S1adj[ site ] , 
			      GAMMAS[ GSRC ] , S1[ site ] ) ;
	}
	corr[ GSRC ][ GSNK ].C[ t ] = (double complex)sum ;
      }
    }
    
    printf("\rContractions done %.f %%",(t+1)/((L0)/100.));fflush(stdout);
  }
  printf("\n");	

  // & do something with the computed correlators
  //#ifdef DEBUG
  debug_mesons( "LL-mesons" , (const struct correlator**)corr ) ;
  //#endif

  // free our correlator measurement
  free_corrs( corr , NSNS , NSNS ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;
  free( S1adj ) ;

  // tell us how long it took, again, probably a long time
  print_time( ) ;

  return SUCCESS ;
}
