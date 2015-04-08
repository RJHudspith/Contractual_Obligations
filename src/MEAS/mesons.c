/**
   @file mesons.c
   @brief meson coputation
 */

#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time() 
#include "io.h"                // read_prop
#include "spinor_ops.h"        // spinor operations
#include "read_propheader.h"   // (re)read the propagator header

// computes meson correlators
int
mesons_diagonal( struct propagator prop ,
		 const char *outfile )
{
  // spinors, S1f is the forward one
  struct spinor *S1 = NULL , *S1f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // data structure for holding the contractions
  struct correlator **wlcorr = NULL , **wwcorr = NULL ;

  // and our spinor
  if( posix_memalign( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( posix_memalign( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    goto free_failure ;
  }

  // allocate the correlators
  wlcorr = allocate_corrs( NSNS , NSNS ) ;
  if( prop.source == WALL ) {
    wwcorr = allocate_corrs( NSNS , NSNS ) ;
  }

  // initially read in a timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto free_failure ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {
    
    // single sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    // above should be non-blocking allowing us to play here
    int GSGK = 0 ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( L0 - 1 ) ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	  //
	}
      }
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < ( NSNS*NSNS ) ; GSGK++ ) {
	
	const int GSRC = GSGK / NSNS ;
	const int GSNK = GSGK % NSNS ;
	
	register double complex sum = 0.0 ;
	
	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S1[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	
	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;
	
	if( prop.source == WALL ) {
	  // correlator computed just out of the summed walls
	  wwcorr[ GSRC ][ GSNK ].C[ t ] =		\
	    meson_contract( GAMMAS[ GSNK ] , SUM1 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
	//
      }
    }

    // to err is human
    if( error_flag == FAILURE ) {
      goto free_failure ;
    }

    // copy S1f into S1
    int i ;
#pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  if( prop.source == WALL ) {
    debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
  }
#endif

  // and write the output files
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)wlcorr ,
		     NSNS , NSNS ) ;

  // free our correlator measurement
  free_corrs( wlcorr , NSNS , NSNS ) ;

  // if we have computed the wall-wall we output it
  if( prop.source == WALL ) {
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)wwcorr ,
		       NSNS , NSNS ) ;
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;
  free( S1f ) ; // free the copy

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

 free_failure :

  // free the correlator
  free_corrs( wlcorr , NSNS , NSNS ) ;

  // free the wall-wall correlator
  if( prop.source == WALL ) {
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ;
  free( S1f ) ; // free the copy

  return FAILURE ;
}

// computes meson correlators
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const char *outfile )
{
  // spinor storage S*f is the forward one being read by the head
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // allocate the gamma basis
  struct gamma *GAMMAS = NULL ;

  // data structure for holding the contractions
  struct correlator **wlcorr = NULL , **wwcorr = NULL ;

  // and our spinors
  if(  posix_memalign( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( posix_memalign( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( posix_memalign( (void**)&S2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( posix_memalign( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
    goto free_failure ;
  }

  // correlator allocations
  wlcorr = allocate_corrs( NSNS , NSNS ) ;
  if( prop1.source == WALL || prop2.source == WALL ) {
    wwcorr = allocate_corrs( NSNS , NSNS ) ;
  }

  // read in the files
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto free_failure ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } else if( prop1.basis == NREL && prop2.basis == CHIRAL ) {
      nrel_rotate_slice( S2 ) ;
    }

    // prop sums
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }

    int error_flag = SUCCESS ;
    int GSGK = 0 ;
    #pragma omp parallel
    {
      // two threads for IO :: TODO, make this more safe
      if( t < ( L0 - 1 ) ) {
        #pragma omp master
	{
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
        #pragma omp single nowait
	{
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }

      // parallelise the furthest out loop
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < ( NSNS * NSNS ) ; GSGK++ ) {
	
	const int GSRC = GSGK / NSNS ;
	const int GSNK = GSGK % NSNS ;

	// accumulate sum
	register double complex sum = 0.0 ;
	
	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S2[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}

	// normal wall-local meson correlator
	wlcorr[ GSRC ][ GSNK ].C[ t ] = sum ;

	if( prop1.source == WALL || prop2.source == WALL ) {
	  // correlator computed just out of the summed walls
	  wwcorr[ GSRC ][ GSNK ].C[ t ] =		\
	    meson_contract( GAMMAS[ GSNK ] , SUM2 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
	//
      }
    }

    // to err is humuan
    if( error_flag == FAILURE ) {
      goto free_failure ;
    }

    // and swap back
    int i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", 100. * (t+1)/(double)(L0) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_mesons( "WL-mesons" , (const struct correlator**)wlcorr ) ;
  if( prop1.source == WALL || prop2.source == WALL ) {
    debug_mesons( "WW-mesons" , (const struct correlator**)wwcorr ) ;
  }
#endif

  // outputs
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)wlcorr ,
		     NSNS , NSNS ) ;
  free_corrs( wlcorr , NSNS , NSNS ) ;

  if( prop1.source == WALL || prop2.source == WALL ) {
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)wwcorr ,
		       NSNS , NSNS ) ;
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S1f ) ;
  free( S2 ) ;
  free( S2f ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

 free_failure :

  // free the walls
  free_corrs( wlcorr , NSNS , NSNS ) ;
  if( prop1.source == WALL || prop2.source == WALL ) {
    free_corrs( wwcorr , NSNS , NSNS ) ;
  }
  
  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;
  free( S1f ) ;
  free( S2 ) ;
  free( S2f ) ;

  return FAILURE ;
}
