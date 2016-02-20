/**
   @file mesons.c
   @brief dispersion relation computation for mesons
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time() 
#include "io.h"                // read_prop
#include "plan_ffts.h"         // ND-1 FFTS
#include "read_propheader.h"   // (re)read the propagator header
#include "setup.h"             // free_ffts() ..
#include "spinor_ops.h"        // sumprop()

// computes flavour-diagonal correlators
int
mesons_diagonal( struct propagator prop ,
		 const struct cut_info CUTINFO ,
		 const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1*stride2 ;

  // spinors, S1f is the forward one
  struct spinor *S1 = NULL , *S1f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // momentum size
  int *NMOM = NULL , *wwNMOM = NULL ;

#ifdef HAVE_FFTW3_H
  // fftw temporary storage
  double complex **in = NULL , **out = NULL ;

  // fftw plans
  fftw_plan *forward = NULL , *backward = NULL ;
#endif

  // allowed momentum list
  struct veclist *list = NULL , *wwlist = NULL;

  // structure containing dispersion relation stuff
  struct mcorr **disp = NULL , **wwdisp = NULL ;

  // loop counters
  size_t i , t ;

  // set a flag we return
  int error_code = SUCCESS ;

  // and our spinor
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 || 
      corr_malloc( (void**)&S1  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

#ifdef HAVE_FFTW3_H

  in  = (double complex**)malloc( flat_dirac * sizeof( double complex* ) ) ;
  out = (double complex**)malloc( flat_dirac * sizeof( double complex* ) ) ;

  forward  = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    in[i]  = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
    out[i] = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
  }

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , flat_dirac , ND-1 ) ;

#endif

  // initialise momentum lists
  init_moms( &NMOM , &wwNMOM , &list , &wwlist , CUTINFO , 
	     prop.source == WALL ? GLU_TRUE : GLU_FALSE ) ;

  // compute the dispersion correlator
  disp = allocate_momcorrs( stride1 , stride2 , NMOM[0] ) ;

  // wall-wall is already zero-projected
  if( prop.source == WALL ) {
    wwdisp = allocate_momcorrs( stride1 , stride2 , wwNMOM[0] ) ;
  }

  // initially read in a timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall-wall sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }
    
    // for multiple time sources
    const size_t tshifted = ( t - prop.origin[ ND-1 ] + LT ) % LT ;

    // master-slave the IO and perform each FFT (if available) in parallel
    size_t GSGK = 0 ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( GAMMAS[ GSNK ] , 
						       GAMMAS[ GAMMA_3 ] ) ;
	// loop spatial hypercube
	size_t site ;
        #ifdef HAVE_FFTW3_H
	for( site = 0 ; site < VOL3 ; site++ ) {
	  in[ GSGK ][ site ] = meson_contract( gt_GSNKdag_gt  , S1[ site ] , 
					       GAMMAS[ GSRC ] , S1[ site ] ,
					       GAMMAS[ GAMMA_5 ] ) ;
	}
	// fft forward
	fftw_execute( forward[ GSGK ] ) ;
	// pack our struct
	size_t p ;
	for( p = 0 ; p < NMOM[0] ; p++ ) {
	  disp[ GSRC ][ GSNK ].mom[ p ].C[ tshifted ] = 
	    out[ GSGK ][ list[ p ].idx ] ;
	}
	#else
	// non - fftw version falls back to the usual zero-momentum sum
	register double complex sum = 0.0 ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( gt_GSNKdag_gt  , S1[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	disp[ GSRC ][ GSNK ].mom[ 0 ].C[ tshifted ] = sum ;
	#endif

	// correlator computed just out of the summed walls
	if( prop.source == WALL ) {
	  wwdisp[ GSRC ][ GSNK ].mom[0].C[ tshifted ] =	\
	    meson_contract( gt_GSNKdag_gt  , SUM1 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
      }
    }

    // to err is human
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy S1f into S1
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)disp , 
		 list , stride1 , stride2 , NMOM ) ;

  // write the wall-wall correlator
  if( prop.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)wwdisp ,
		   wwlist , stride1 , stride2 , wwNMOM ) ;
  }

 memfree :

  // free correlators and momentum list
  if( NMOM != NULL ) {
    free_momcorrs( disp , stride1 , stride2 , NMOM[0] ) ;
    if( prop.source == WALL ) {
      free_momcorrs( wwdisp , stride1 , stride2 , wwNMOM[0] ) ;
    }
  }

  // free our ffts
  free_ffts( in , out , forward , backward , flat_dirac ) ;

  // free momenta lists
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}
