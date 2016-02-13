/**
   @file tetra_degen.c
   @brief tetraquark contraction code for degenerate light quarks
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

// tetraquark candidate L1 L1 \bar{H} \bar{H} so two propagators
int
tetraquark_degen( struct propagator prop1 ,
		  struct propagator prop2 ,
		  struct cut_info CUTINFO ,
		  const char *outfile )
{
  // flat dirac indices
  const size_t flat_dirac = TETRA_NOPS * ( B_CHANNELS ) ; // is gamma_i

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
  
  // loop counter
  size_t i , t ;

  // error code
  int error_code = SUCCESS ;

  // allocations
  if( corr_malloc( (void**)&S1  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( prop1.basis == NREL || prop2.basis == NREL ) { 
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      error_code = FAILURE ; goto memfree ;
    }
  } else {
    if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      error_code = FAILURE ; goto memfree ;
    }
  }

  NMOM = malloc( sizeof( int ) ) ;
  wwNMOM = malloc( sizeof( int ) ) ;

  in = malloc( flat_dirac * sizeof( double complex* ) ) ;
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

  list = (struct veclist*)compute_veclist( NMOM , CUTINFO , 
					   ND-1 , GLU_FALSE ) ;

#else

  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

#endif

  // allocate the zero veclist for our writing out purposes
  if( prop1.source == WALL || prop2.source == WALL ) {
    wwlist = (struct veclist*)zero_veclist( wwNMOM , ND-1 , GLU_FALSE ) ;
  }

  // Output with TETRA_OPS channels and B_CHANNELS 
  tetra_corr = allocate_momcorrs( TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL || prop2.source == WALL ) {
    tetra_corrWW = allocate_momcorrs( TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE || 
      read_prop( prop2 , S2 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

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
    struct spinor SUM1 , SUM2 , SUMbwdH ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
      full_adj( &SUMbwdH , SUM2 , GAMMAS[ GAMMA_5 ] ) ;
    }

    // assumes all sources are at the same origin, checked in wrap_tetras
    const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 

    // strange memory access pattern threads better than what was here before
    size_t site ;
    #pragma omp parallel
    {
      // read on the master and one slave
      if( t < LT-1 ) {
        #pragma omp master
	{
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
        #pragma omp single nowait
	{
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// precompute backward bottom propagator using 
	// gamma_5 hermiticity
	struct spinor bwdH ;
	full_adj( &bwdH , S2[ site ] , GAMMAS[ GAMMA_5 ] ) ;

	// diquark-diquark tetra
	double complex result[ TETRA_NOPS ] ;

	// loop gamma source
	size_t GSRC , op ;
	for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
	  // perform contraction, result in result
	  tetras( result , S1[ site ] , S1[ site ] , bwdH , GAMMAS , GSRC , 
		  GLU_TRUE ) ;
	  // put contractions into flattend array for FFT
	  for( op = 0 ; op < TETRA_NOPS ; op++ ) {
	    in[ op + TETRA_NOPS * GSRC ][ site ] = result[ op ] ;
	  }
	}
      }
      // wall-wall contractions
      if( prop1.source == WALL || prop2.source == WALL ) {
	double complex result[ TETRA_NOPS ] ;
	size_t GSRC , op ;
	for( GSRC = 0 ; GSRC < B_CHANNELS ; GSRC++ ) {
	  // perform contraction, result in result
	  tetras( result , SUM1 , SUM1 , SUMbwdH , GAMMAS , GSRC , GLU_TRUE ) ;
	  // put contractions into final correlator object
	  for( op = 0 ; op < TETRA_NOPS ; op++ ) {
	    tetra_corrWW[ op ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[ op ] ;
	  }
	}
      }
      // end of walls
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
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf( "\r[TETRA] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)tetra_corr ,
		 list , TETRA_NOPS , B_CHANNELS , NMOM ) ;

  // if we have walls we use them
  if( prop1.source == WALL || prop2.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)tetra_corrWW ,
		   wwlist , TETRA_NOPS , B_CHANNELS , wwNMOM ) ;
  }

  // failure sink
 memfree :

  // free our correlators
  if( NMOM != NULL ) {
    free_momcorrs( tetra_corr , TETRA_NOPS , B_CHANNELS , NMOM[0] ) ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      free_momcorrs( tetra_corrWW , TETRA_NOPS , B_CHANNELS , wwNMOM[0] ) ;
    }
  }

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

  // free spinors
  free( S1f ) ; free( S1 ) ; 
  free( S2f ) ; free( S2 ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}

