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

// tetraquark candidate L1 L2 \bar{H} \bar{H} so three propagators
int
tetraquark( struct propagator prop1 ,
	    struct propagator prop2 ,
	    struct propagator prop3 ,
	    struct cut_info CUTINFO ,
	    const char *outfile )
{
  // flat dirac indices
  const size_t flat_dirac = TETRA_NOPS * ( B_CHANNELS ) ; // is gamma_i

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL , 
    *S2 = NULL , *S2f = NULL , *S3 = NULL , *S3f = NULL ;

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
  if( corr_malloc( (void**)&S3  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S3f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( prop1.basis == NREL || prop2.basis == NREL || prop3.basis == NREL ) { 
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
  if( prop1.source == WALL || prop2.source == WALL || prop3.source == WALL ) {
    wwlist = (struct veclist*)zero_veclist( wwNMOM , B_CHANNELS , GLU_FALSE ) ;
  }

  // Define our output correlators, with 2 channels and ND-1 components
  tetra_corr = allocate_momcorrs( TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL || prop2.source == WALL || prop3.source == WALL ) {
    tetra_corrWW = allocate_momcorrs( TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE || 
      read_prop( prop2 , S2 ) == FAILURE || 
      read_prop( prop3 , S3 )  == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral hadrons we switch chiral to nrel
    if( prop1.basis == CHIRAL && ( prop2.basis == NREL || prop3.basis == NREL ) ) {
      nrel_rotate_slice( S1 ) ;
    } 
    if( prop2.basis == CHIRAL && ( prop1.basis == NREL || prop3.basis == NREL ) ) {
      nrel_rotate_slice( S2 ) ;
    } 
    if( prop3.basis == CHIRAL && ( prop1.basis == NREL || prop2.basis == NREL ) ) {
      nrel_rotate_slice( S3 ) ;
    }

    // compute wall sum
    struct spinor SUM1 , SUM2 , SUM3 , SUMbwdH ;
    if( prop1.source == WALL || prop2.source == WALL || prop3.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
      sumprop( &SUM3 , S3 ) ;
      full_adj( &SUMbwdH , SUM2 , GAMMAS[ GAMMA_5 ] ) ;
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
        #pragma omp single nowait
	{
	  if( read_prop( prop3 , S3f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// precompute backward bottom propagator
	struct spinor bwdH ;
	full_adj( &bwdH , S3[ site ] , GAMMAS[ GAMMA_5 ] ) ;

	// diquark-diquark tetra
	double complex result[ TETRA_NOPS - 1 ] ;

	// loop gamma source
	size_t GSRC ;
	for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
	  // poke in the tetra
	  if( diquark_diquark( result , S1[ site ] , S2[ site ] , bwdH , 
			       GAMMAS , GSRC ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	  // for each C_\gamma_i we put out all NOPS operators
	  in[ 0 + TETRA_NOPS*GSRC ][ site ] = result[0] ; // O1 O1^\dagger
	  in[ 1 + TETRA_NOPS*GSRC ][ site ] = result[1] ; // O1 O2^\dagger
	  in[ 2 + TETRA_NOPS*GSRC ][ site ] = result[2] ; // O2 O1^\dagger
	  in[ 3 + TETRA_NOPS*GSRC ][ site ] = result[3] ; // O2 O2^\dagger
	  // poke in the dimeson part :: O3 O3^\dagger
	  in[ 4 + TETRA_NOPS*GSRC ][ site ] = 
	    2 * dimeson( S1[ site ] , S2[ site ] , bwdH , GAMMAS , GSRC ) ;
	}
      }
      // wall-wall contractions
      if( prop1.source == WALL || prop2.source == WALL ) {
	double complex result[ TETRA_NOPS - 1 ] ;
	size_t GSRC ;
	for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
	  // compute diquark-diquark tetra
	  if( diquark_diquark( result , SUM1 , SUM2 , SUMbwdH , 
			       GAMMAS , GSRC ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	  // diquark-diquark
	  tetra_corrWW[ 0 + TETRA_NOPS*GSRC ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[0] ;
	  tetra_corrWW[ 1 + TETRA_NOPS*GSRC ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[1] ;
	  tetra_corrWW[ 2 + TETRA_NOPS*GSRC ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[2] ;
	  tetra_corrWW[ 3 + TETRA_NOPS*GSRC ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[3] ;
	  // meson-meson
	  tetra_corrWW[ 4 + TETRA_NOPS*GSRC ][ GSRC ].mom[ 0 ].C[ tshifted ] = 
	    dimeson( SUM1 , SUM2 , SUMbwdH , GAMMAS , GSRC ) +
	    dimeson( SUM2 , SUM1 , SUMbwdH , GAMMAS , GSRC ) ;
	}
	///
      }
    }

    // momentum projection
    size_t GSRC ;
    #pragma omp parallel for private(GSRC) schedule(dynamic)
    for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
      size_t op , p ;
      for( op = 0 ; op < TETRA_NOPS ; op++ ) {
        #ifdef HAVE_FFTW3_H
	fftw_execute( forward[ op + TETRA_NOPS * GSRC ] ) ;
	for( p = 0 ; p < NMOM[0] ; p++ ) {
	  tetra_corr[ op ][ GSRC ].mom[ p ].C[ tshifted ] =
	    out[ op + TETRA_NOPS*GSRC ][ list[ p ].idx ] ;
	}
        #else
	register double complex sum = 0.0 ;
	for( p = 0 ; p < LCU ; p++ ) {
	  sum += in[ op + TETRA_NOPS*GSRC ][ p ] ;
	}
	tetra_corr[ op ][ GSRC ].mom[ 0 ].C[ tshifted ] = sum ;
        #endif
      }
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
      memcpy( &S3[i] , &S3f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[TETRA] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)tetra_corr ,
		 list , TETRA_NOPS , B_CHANNELS , NMOM ) ;
  free_momcorrs( tetra_corr , TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;

  // if we have walls we use them
  if( prop1.source == WALL || prop2.source == WALL || prop3.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)tetra_corrWW ,
		   wwlist , TETRA_NOPS , B_CHANNELS , wwNMOM ) ;
    free_momcorrs( tetra_corrWW , TETRA_NOPS , B_CHANNELS , wwNMOM[0] ) ;
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
  free( S3 ) ; free( S3f ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;
  rewind( prop3.file ) ; read_propheader( &prop3 ) ;

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
    free_momcorrs( tetra_corr , TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;
    if( prop1.source == WALL || prop2.source == WALL || prop3.source == WALL ) {
      free_momcorrs( tetra_corrWW , TETRA_NOPS , B_CHANNELS , wwNMOM[0] ) ;
    }
  }

  // free spinors
  free( S1f ) ; free( S1 ) ; 
  free( S2f ) ; free( S2 ) ;
  free( S3f ) ; free( S3 ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}

