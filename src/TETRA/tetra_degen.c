/**
   @file tetraquark.c
   @brief tetraquark contraction code
*/

#include "common.h"

#include "basis_conversions.h"  // nrel_rotate_slice()
#include "contractions.h"       // gamma_mul_lr()
#include "correlators.h"        // allocate_corrs() && free_corrs()
#include "cut_routines.h"       // veclist
#include "gammas.h"             // make_gammas() && gamma_mmul*
#include "GLU_timer.h"          // print_time()
#include "io.h"                 // for read_prop()
#include "plan_ffts.h"          // create_plans_DFT() 
#include "read_propheader.h"    // for read_propheader()
#include "spinor_ops.h"         // sumprop()
#include "tetra_contractions.h" // diquark_diquark()

// tetraquark candidate L1 L2 \bar{H} \bar{H} but we can have the 
// case where L1 == L2 say for ud components
int
tetraquark_degen( struct propagator prop1 ,
		  struct propagator prop2 ,
		  struct cut_info CUTINFO ,
		  const char *outfile )
{
  // flat dirac indices
  const int flat_dirac = 2 * ( ND-1 ) ; // is gamma_i

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // momentum stuff
  int *NMOM = NULL , *wwNMOM = NULL ;
  struct veclist *list = NULL , *wwlist = NULL ;

  // fftw temporaries
  double complex **in = NULL , **out = NULL ;

  // correlators
  struct mcorr **tetra_corr = NULL , **tetra_corrWW = NULL ;

#ifdef HAVE_FFTW3_H
  fftw_plan *forward = NULL , *backward = NULL ;
#else
  int *forward = NULL , *backward = NULL ;
#endif

  // allocations
  if( corr_malloc( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S2  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( prop1.basis == NREL || prop2.basis == NREL ) { 
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      goto FREE_FAIL ;
    }
  } else {
    if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      goto FREE_FAIL ;
    }
  }

  NMOM = malloc( sizeof( int ) ) ;
  wwNMOM = malloc( sizeof( int ) ) ;

  in = malloc( flat_dirac * sizeof( double complex* ) ) ;
  size_t i ;
  for( i = 0 ; i < flat_dirac ; i++ ) {
    in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H

  out = malloc( flat_dirac * sizeof( double complex* ) ) ;
  for( i = 0 ; i < flat_dirac ; i++ ) {
    out[ i ] = malloc( LCU * sizeof( double complex ) ) ;
  }

  forward  = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ;

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , flat_dirac , ND-1 ) ;

  list = (struct veclist*)compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

#else

  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

#endif

  // allocate the zero veclist for our writing out purposes
  if( prop1.source == WALL || prop2.source == WALL ) {
    wwlist = (struct veclist*)zero_veclist( wwNMOM , ND-1 , GLU_FALSE ) ;
  }

  // Define our output correlators, with 2 channels and ND-1 components
  tetra_corr = allocate_momcorrs( 2 , ND-1 , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL || prop2.source == WALL ) {
    tetra_corrWW = allocate_momcorrs( 2 , ND-1 , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE || 
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral hadrons we switch chiral to nrel
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } 
    if( prop2.basis == CHIRAL && prop1.basis == NREL ) {
      nrel_rotate_slice( S2 ) ;
    } 

    // compute wall sum
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }

    // assumes all sources are at the same origin, checked in wrap_tetras
    const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 

    // strange memory access pattern threads better than what was here before
    size_t site ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      // read on the master and one slave
      if( t < LT-1 ) {
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
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// precompute backward bottom propagator
	struct spinor bwdH ;
	full_adj( &bwdH , S2[ site ] , GAMMAS[ GAMMA_5 ] ) ;

	// loop gamma source
	size_t GSRC ;
	for( GSRC = 0 ; GSRC < ND-1 ; GSRC++ ) {
	  // poke in the tetra
	  in[ GSRC ][ site ] = diquark_diquark( S1[ site ] , S1[ site ] , bwdH , 
						GAMMAS , GSRC ) ;
	  // poke in the dimesons
	  in[ GSRC + ND-1 ][ site ] = 
	    2.0 * dimeson( S1[ site ] , S1[ site ] , bwdH , GAMMAS , GSRC ) ;
	}
      }
      // wall-wall contractions
      if( prop1.source == WALL || prop2.source == WALL ) {
	struct spinor bwdH ;
	full_adj( &bwdH , SUM2 , GAMMAS[ GAMMA_5 ] ) ;
	size_t GSRC ;
	for( GSRC = 0 ; GSRC < ND-1 ; GSRC++ ) {
	  // diquark-diquark
	  tetra_corrWW[ 0 ][ GSRC ].mom[ 0 ].C[ tshifted ] = 
	    diquark_diquark( SUM1 , SUM1 , bwdH , GAMMAS , GSRC ) ;
	  // meson-meson
	  tetra_corrWW[ 1 ][ GSRC ].mom[ 0 ].C[ tshifted ] = 
	    2.0 * dimeson( SUM1 , SUM1 , bwdH , GAMMAS , GSRC ) ;
	}
	///
      }
    }

    // momentum projection
    size_t GSRC , p ;
    #pragma omp parallel for private(GSRC) schedule(dynamic)
    for( GSRC = 0 ; GSRC < ND-1 ; GSRC++ ) {
      #ifdef HAVE_FFTW3_H
      fftw_execute( forward[ GSRC ] ) ;
      fftw_execute( forward[ GSRC + ND-1 ] ) ;
      for( p = 0 ; p < NMOM[0] ; p++ ) {
	tetra_corr[ 0 ][ GSRC ].mom[ p ].C[ tshifted ] =
	  out[ GSRC ][ list[ p ].idx ] ;
	tetra_corr[ 1 ][ GSRC ].mom[ p ].C[ tshifted ] =
	  out[ GSRC + ND-1 ][ list[ p ].idx ] ;
      }
      #else
      register double complex sum1 , sum2 ;
      for( p = 0 ; p < LCU ; p++ ) {
	sum1 += in[ GSRC ][ p ] ;
	sum2 += in[ GSRC + ND - 1 ][ p ] ;
      }
      tetra_corr[ 0 ][ GSRC ].mom[ 0 ].C[ tshifted ] = sum1 ;
      tetra_corr[ 1 ][ GSRC ].mom[ 0 ].C[ tshifted ] = sum2 ;
      #endif
    }

    // if we error we leave
    if( error_flag == FAILURE ) {
      goto FREE_FAIL ;
    }

    // copy over the propagators
    size_t i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[TETRA] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)tetra_corr ,
		 list , 2 , ND-1 , NMOM ) ;
  free_momcorrs( tetra_corr , 2 , ND-1 , NMOM[0] ) ;

  // if we have walls we use them
  if( prop1.source == WALL || prop2.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)tetra_corrWW ,
		   wwlist , 2 , ND-1 , wwNMOM ) ;
    free_momcorrs( tetra_corrWW , 2 , ND-1 , wwNMOM[0] ) ;
  }

  // free the "in" allocation
  for( i = 0 ; i < flat_dirac ; i++ ) {
    free( in[ i ] ) ;
  }
  free( in ) ;

#ifdef HAVE_FFTW3_H
  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < flat_dirac ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; 
  fftw_cleanup( ) ; 
#endif

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free stuff
  free( S1 ) ; free( S1f ) ;
  free( S2 ) ; free( S2f ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

  // failure sink
 FREE_FAIL :

  // free the "in" allocation
  if( in != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      free( in[ i ] ) ;
    }
  }
  free( in ) ;

#ifdef HAVE_FFTW3_H
  // free fftw stuff
  if( out != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      free( out[ i ] ) ;
    }
  }
  if( forward != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_destroy_plan( forward[ i ] ) ;
    }
  }
  if( backward != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_destroy_plan( forward[ i ] ) ;
    }
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; 
  fftw_cleanup( ) ; 
#endif

  // free our correlators
  if( NMOM != NULL ) {
    free_momcorrs( tetra_corr , 2 , ND-1 , NMOM[0] ) ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      free_momcorrs( tetra_corrWW , 2 , ND-1 , wwNMOM[0] ) ;
    }
  }

  // free spinors
  free( S1f ) ; free( S1 ) ; 
  free( S2f ) ; free( S2 ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}

