/**
   @file baryons.c
   @brief baryon contraction code
 */

#include "common.h"

#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // meson contract or whatever
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas()
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()

int
baryons_diagonal( struct propagator prop ,
		  const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // and our spinor
  struct spinor *S1 ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free( S1 ) ; free( GAMMAS ) ;
    printf( "[MESONS] memalign failure \n" ) ;
    return FAILURE ;
  }

  struct correlator **corr = allocate_corrs( NSNS , NSNS ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop , S1 ) == FAILURE ) {
      free_corrs( corr , NSNS , NSNS ) ; 
      free( GAMMAS ) ; free( S1 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop.basis == CHIRAL && prop.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } 

    //
    // Contractions go here
    // 

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {

      }
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)corr ,
		     NSNS , NSNS ) ;

  free_corrs( corr , NSNS , NSNS ) ;

  free( S1 ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}

