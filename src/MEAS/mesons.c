/**
   @file mesons.c
   @brief dispersion relation computation for mesons
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "cut_routines.h"      // momentum list and stuff
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time() 
#include "io.h"                // read_prop
#include "plan_ffts.h"         // ND-1 FFTS
#include "read_propheader.h"   // (re)read the propagator header
#include "spinor_ops.h"        // sumprop()

// computes flavour-diagonal correlators
int
mesons_diagonal( struct propagator prop ,
		 const struct cut_info CUTINFO ,
		 const char *outfile )
{
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

  int i ;

  // and our spinor
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    goto free_failure ;
  }

  // compute the momentum list for the specified cut
  NMOM = (int*)malloc( sizeof( int ) ) ;

#ifdef HAVE_FFTW3_H

  in  = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;
  out = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;

  forward  = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    in[i]  = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
    out[i] = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
  }

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , ( NSNS*NSNS ) , ND-1 ) ;

  list = (struct veclist*)compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

#else

  // create a ( 0 , 0 , 0 ) vector list
  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

#endif

  // compute the dispersion correlator
  disp = allocate_momcorrs( NSNS , NSNS , NMOM[0] ) ;

  // wall-wall is already zero-projected
  if( prop.source == WALL ) {
    wwNMOM = (int*)malloc( sizeof( int ) ) ;
    wwlist = (struct veclist*)zero_veclist( wwNMOM , ND-1 , GLU_FALSE ) ;
    wwdisp = allocate_momcorrs( NSNS , NSNS , wwNMOM[0] ) ;
  }

  // initially read in a timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto free_failure ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // compute wall-wall sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    // master-slave the IO and perform each FFT (if available) in parallel
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
	// loop spatial hypercube
	int site ;
        #ifdef HAVE_FFTW3_H
	for( site = 0 ; site < VOL3 ; site++ ) {
	  in[ GSGK ][ site ] = meson_contract( GAMMAS[ GSNK ] , S1[ site ] , 
					       GAMMAS[ GSRC ] , S1[ site ] ,
					       GAMMAS[ GAMMA_5 ] ) ;
	}
	// fft forward
	fftw_execute( forward[ GSGK ] ) ;
	// pack our struct
	int p ;
	for( p = 0 ; p < NMOM[0] ; p++ ) {
	  disp[ GSRC ][ GSNK ].mom[ p ].C[ t ] = out[ GSGK ][ list[ p ].idx ] ;
	}
	#else
	// non - fftw version falls back to the usual zero-momentum sum
	register double complex sum = 0.0 ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  sum += meson_contract( GAMMAS[ GSNK ] , S1[ site ] , 
				 GAMMAS[ GSRC ] , S1[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ;
	}
	disp[ GSRC ][ GSNK ].mom[ 0 ].C[ t ] = sum ;
	#endif

	// correlator computed just out of the summed walls
	if( prop.source == WALL ) {
	  wwdisp[ GSRC ][ GSNK ].mom[0].C[ t ] =	\
	    meson_contract( GAMMAS[ GSNK ] , SUM1 , 
			    GAMMAS[ GSRC ] , SUM1 ,
			    GAMMAS[ GAMMA_5 ] ) ;
	}
      }
    }

    // to err is human
    if( error_flag == FAILURE ) {
      goto free_failure ;
    }

    // copy S1f into S1
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[MESONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)disp , list , 
		 NSNS , NSNS , NMOM ) ;

  // free momenta for the wall-local
  free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;

  // write the wall-wall correlator
  if( prop.source == WALL ) {
    char outstr[ 256 ] ;
    sprintf( outstr , "%s.ww" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)wwdisp ,
		   wwlist , NSNS , NSNS , wwNMOM ) ;
    free_momcorrs( wwdisp , NSNS , NSNS , wwNMOM[0] ) ;
  }
  
#ifdef HAVE_FFTW3_H
  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
    fftw_free( in[i] ) ;
  }
  free( forward ) ; 
  free( backward ) ; 
  fftw_free( out ) ; 
  fftw_free( in ) ; 
  fftw_cleanup( ) ; 
#endif

  // free number of possible momenta and the momentum list
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

  return SUCCESS ;

 free_failure :

#ifdef HAVE_FFTW3_H
  if( in != NULL ) {
    for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
      free( in[i] ) ;
    }
  }
  if( out != NULL ) {
    for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
      free( out[i] ) ;
    }
  }
  if( forward != NULL ) {
    for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
      fftw_destroy_plan( forward[i] ) ;
    }
  }
  if( backward != NULL ) {
    for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
      fftw_destroy_plan( backward[i] ) ;
    }
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; fftw_free( in ) ; 
  fftw_cleanup( ) ; 
#endif

  // free correlators and momentum list
  if( NMOM != NULL ) {
    free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;
    if( prop.source == WALL ) {
      free_momcorrs( wwdisp , NSNS , NSNS , wwNMOM[0] ) ;
    }
  }

  // free momenta lists
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor(s)
  free( S1 ) ; free( S1f ) ;

  return FAILURE ;
}
