/**
   @file mesons_offdiag.c
   @brief flavour off-diagonal meson contractions

   split from mesons for clarity
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

// computes flavour-off diagonal correlators with ND-1 momentum projection
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const struct cut_info CUTINFO ,
		    const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1*stride2 ;

  // spinors, S1f is the forward one
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // number of allowed momenta
  int *NMOM = NULL , *wwNMOM = NULL ;

#ifdef HAVE_FFTW3_H
  // fftw temporary storage
  double complex **in = NULL , **out = NULL ;

  // fftw plans
  fftw_plan *forward = NULL , *backward = NULL ;
#endif

  // set momenta and the momenta list
  struct veclist *list = NULL , *wwlist = NULL ;

  // structure containing dispersion relation stuff
  struct mcorr **disp = NULL , **wwdisp = NULL ;

  // error code
  int error_code = SUCCESS ;

  // loop counters
  size_t i , t ;

  // and our spinor
  if( corr_malloc( (void**)&S1  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
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
	     prop1.source == WALL || prop2.source == WALL ? \
	     GLU_TRUE : GLU_FALSE ) ;

  // compute the dispersion correlator
  disp = allocate_momcorrs( stride1 , stride2 , NMOM[0] ) ;

  // and the walls
  if( prop1.source == WALL || prop2.source == WALL ) {
    wwdisp = allocate_momcorrs( stride1 , stride2 , wwNMOM[0] ) ;
  }

  // read in the files
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop1.basis == CHIRAL && ( prop2.basis == NREL ) ) {
      nrel_rotate_slice( S1 ) ;
    } else if( ( prop1.basis == NREL ) && prop2.basis == CHIRAL ) {
      nrel_rotate_slice( S2 ) ;
    }

    // prop sums
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL || prop2.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }
    
    // support for multiple time sources
    const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ;

    // master-slave the IO and perform each FFT in parallel
    size_t GSGK ;
    #pragma omp parallel
    {
     // two threads for IO
      if( t < ( LT - 1 ) ) {
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
	  in[ GSGK ][ site ] = meson_contract( gt_GSNKdag_gt  , S2[ site ] ,
					       GAMMAS[ GSRC ] , S1[ site ] ,
					       GAMMAS[ GAMMA_5 ] ) ;
	}
	// fft forward is e^{ip.x}
	fftw_execute( forward[ GSGK ] ) ;
	// pack our struct
	size_t p ;
	for( p = 0 ; p < NMOM[0] ; p++ ) {
	  disp[ GSRC ][ GSNK ].mom[ p ].C[ tshifted ] =
	    out[ GSGK ][ list[ p ].idx ] ;
	}
        #else
	register double complex sum = 0.0 ;
	// for the non-fftw'd version we fall back on the zero-mom projection
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( gt_GSNKdag_gt  , S2[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	disp[ GSRC ][ GSNK ].mom[ 0 ].C[ tshifted ] = sum ;
	#endif

	// and contract the walls
	if( prop1.source == WALL || prop2.source == WALL ) {
	  wwdisp[ GSRC ][ GSNK ].mom[ 0 ].C[ tshifted ] = \
	    meson_contract( gt_GSNKdag_gt  , SUM2 ,
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
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)disp , list , 
		 stride1 , stride2 , NMOM ) ;

  // and do the walls as they are basically free
  if( prop1.source == WALL || prop2.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)wwdisp ,
		   wwlist , stride1 , stride2 , wwNMOM ) ;
  }

  // memory freeing part
 memfree :

  // free correlators and momentum list
  if( NMOM != NULL ) {
    free_momcorrs( disp , stride1 , stride2 , NMOM[0] ) ;
    if( prop1.source == WALL ) {
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

  // free our spinors
  free( S1 ) ;  free( S1f ) ;
  free( S2 ) ;  free( S2f ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}
