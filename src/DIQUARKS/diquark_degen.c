/**
   @file diquark.c
   @brief diquark-diquark contraction
 */
#include "common.h"

#include "basis_conversions.h"    // nrel_rotate_slice()
#include "correlators.h"          // allocate_corrs() && free_corrs()
#include "diquark_contraction.h"  // diquark()
#include "gammas.h"               // make_gammas() && gamma_mmul*
#include "GLU_timer.h"            // print_time()
#include "io.h"                   // for read_prop()
#include "plan_ffts.h"            // create_plans_DFT() 
#include "read_propheader.h"      // for read_propheader()
#include "setup.h"                // compute_correlator() ..
#include "spinor_ops.h"           // sumprop()

// degenerate diquarks
int
diquark_degen( struct propagator prop1 ,
	       const struct cut_info CUTINFO , 
	       const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1*stride2 ;

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL ;

  // momentum stuff
  int *NMOM = NULL , *wwNMOM = NULL ;
  struct veclist *list = NULL , *wwlist = NULL ;

  // fftw temporaries
  double complex **in = NULL , **out = NULL ;

  // correlators
  struct mcorr **diquark_corr = NULL , **diquark_corrWW = NULL ;

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
      corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop1.basis ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // create a look up table for the gammas, these are square 
  // so we only do the first NSNS
  struct gamma *GAM1 = malloc( stride1 * sizeof( struct gamma ) ) ;
  struct gamma *GAM2 = malloc( stride2 * sizeof( struct gamma ) ) ; 
  for( i = 0 ; i < stride1 ; i++ ) {
    GAM1[ i ]= CGmu( GAMMAS[ i ] , GAMMAS ) ;
  }
  for( i = 0 ; i < stride2 ; i++ ) {
    GAM2[ i ]= gt_Gdag_gt( CGmu( GAMMAS[ i ] , GAMMAS ) , 
			   GAMMAS[ GAMMA_3 ] ) ;
  }

  // allocate result array
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

#endif

  // initialise momentum lists
  init_moms( &NMOM , &wwNMOM , &list , &wwlist , CUTINFO , 
	     prop1.source == WALL ? GLU_TRUE : GLU_FALSE ) ;

  // Output with stride1*stride2 size correlators
  diquark_corr = allocate_momcorrs( stride1 , stride2 , NMOM[0] ) ;

  // allocate the wall correlators if we are using wall source propagators
  if( prop1.source == WALL ) {
    diquark_corrWW = allocate_momcorrs( stride1 , stride2 , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall sum
    struct spinor SUM1 ;
    if( prop1.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
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
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// loop gamma source
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1*stride2 ; GSGK++ ) {
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  // perform contraction, result in result
	  in[ GSNK + stride2*GSRC ][ site ] = 
	    -2 * diquark( S1[ site ] , S1[ site ] , 
			  GAM1[ GSRC ] , GAM2[ GSNK ] ) ;
	}
      }
      // wall-wall contractions
      if( prop1.source == WALL ) {
	// loop gamma source
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1*stride2 ; GSGK++ ) {
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  diquark_corrWW[ GSRC ][ GSNK ].mom[0].C[ tshifted ] = 
	    diquark( SUM1 , SUM1 , GAM1[ GSRC ] , GAM2[ GSNK ] ) ;
	}
      }
      // end of walls
    }

    // function computes the correlator, fftw-ing if available
    compute_correlator( diquark_corr , 
			(const double complex**)in , 
			(const double complex**)out , 
			list , NMOM , forward , stride1 , stride2 , 
			tshifted ) ;

    // if we error we leave
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf( "\r[DIQUARK] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)diquark_corr ,
		 list , stride1 , stride2 , NMOM ) ;

  // if we have walls we use them
  if( prop1.source == WALL ) {
    printf( "Walls not done yet \n" ) ;
  }

  // failure sink
 memfree :

  // free our correlators
  if( NMOM != NULL ) {
    free_momcorrs( diquark_corr , stride1 , stride2 , NMOM[0] ) ;
    if( prop1.source == WALL ) {
      free_momcorrs( diquark_corrWW , stride1 , stride2 , wwNMOM[0] ) ;
    }
  }

  // free our ffts
  free_ffts( in , out , forward , backward , flat_dirac ) ;

  // free spinors
  free( S1f ) ; free( S1 ) ; 

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ; free( GAM1 ) ; free( GAM2 ) ;

  // rewind file and read header again, this should be in the 
  // wrapper no?
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}
