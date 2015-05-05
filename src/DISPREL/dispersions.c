/**
   @file dispersions.c
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

// computes flavour-diagonal correlators
int
dispersions_diagonal( struct propagator prop ,
		      const struct cut_info CUTINFO ,
		      const char *outfile )
{
#ifndef HAVE_FFTW3_H
  printf( "[DISPREL] routines are only available when FFTW is linked! \n" ) ;
  return FAILURE ;
#else
  // spinors, S1f is the forward one
  struct spinor *S1 = NULL , *S1f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // fftw temporary storage
  double complex **in = NULL , **out = NULL ;

  // fftw plans
  fftw_plan *forward = NULL , *backward = NULL ;

  // set momenta and the momenta list
  int *NMOM = NULL ;

  // structure containing dispersion relation stuff
  struct mcorr **disp = NULL ;

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

  in  = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;
  out = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;

  forward  = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    in[i]  = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
    out[i] = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
  }

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , ( NSNS*NSNS ) , ND-1 ) ;

  // compute the momentum list for the specified cut
  NMOM = (int*)malloc( sizeof( int ) ) ;
  const struct veclist *list = compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

  // compute the dispersion correlator
  disp = allocate_momcorrs( NSNS , NSNS , NMOM[0] ) ;

  // initially read in a timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto free_failure ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // master-slave the IO and perform each FFT in parallel
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
    printf("\r[DISPREL] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)disp , list , 
		 NSNS , NSNS , NMOM ) ;

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

  // free momenta
  free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;

  free( NMOM ) ;
  free( (void*)list ) ;

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

  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
    fftw_free( in[i] ) ;
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; fftw_free( in ) ; 
  fftw_cleanup( ) ; 

  // free correlators and momentum list
  free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;

  free( NMOM ) ;
  free( (void*)list ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinor
  free( S1 ) ; free( S1f ) ; 
  return FAILURE ;
#endif
}

// computes flavour-off diagonal correlators with ND-1 momentum projection
int
dispersions_offdiagonal( struct propagator prop1 ,
			 struct propagator prop2 ,
			 const struct cut_info CUTINFO ,
			 const char *outfile )
{
#ifndef HAVE_FFTW3_H
  printf( "[DISPREL] routines are only available when FFTW is linked! \n" ) ;
  return FAILURE ;
#else
  // spinors, S1f is the forward one
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = NULL ;

  // fftw temporary storage
  double complex **in = NULL , **out = NULL ;

  // fftw plans
  fftw_plan *forward = NULL , *backward = NULL ;

  // set momenta and the momenta list
  int *NMOM = NULL ;

  // structure containing dispersion relation stuff
  struct mcorr **disp = NULL ;

  // and our spinor
  if( posix_memalign( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
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
  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , ( prop1.basis == NREL || prop2.basis == NREL ) ? \
		   NREL : prop1.basis ) == FAILURE ) {
    goto free_failure ;
  }

  in  = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;
  out = (double complex**)malloc( ( NSNS*NSNS ) * sizeof( double complex* ) ) ;

  forward  = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( ( NSNS*NSNS ) * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    in[i]  = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
    out[i] = ( double complex* )fftw_malloc( LCU * sizeof( double complex ) ) ; 
  }

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , ( NSNS*NSNS ) , ND-1 ) ;

  // compute the momentum list for the specified cut
  NMOM = (int*)malloc( sizeof( int ) ) ;
  const struct veclist *list = compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

  // compute the dispersion correlator
  disp = allocate_momcorrs( NSNS , NSNS , NMOM[0] ) ;

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

    // master-slave the IO and perform each FFT in parallel
    int GSGK = 0 ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
     // two threads for IO
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
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < ( NSNS*NSNS ) ; GSGK++ ) {
	const int GSRC = GSGK / NSNS ;
	const int GSNK = GSGK % NSNS ;
	// loop spatial hypercube
	int site ;
	for( site = 0 ; site < VOL3 ; site++ ) {
	  in[ GSGK ][ site ] = meson_contract( GAMMAS[ GSNK ] , S2[ site ] , 
					       GAMMAS[ GSRC ] , S1[ site ] ,
					       GAMMAS[ GAMMA_5 ] ) ;
	}
	// fft forward is e^{ip.x}
	fftw_execute( forward[ GSGK ] ) ;

	// pack our struct
	int p ;
	for( p = 0 ; p < NMOM[0] ; p++ ) {
	  disp[ GSRC ][ GSNK ].mom[ p ].C[ t ] = out[ GSGK ][ list[ p ].idx ] ;
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
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[DISPREL] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)disp , list , 
		 NSNS , NSNS , NMOM ) ;

  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
    fftw_free( in[i] ) ;
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; fftw_free( in ) ; 
  fftw_cleanup( ) ; 

  // free momenta
  free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;

  free( NMOM ) ;
  free( (void*)list ) ;

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

  return SUCCESS ;

 free_failure :

  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( NSNS*NSNS ) ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
    fftw_free( in[i] ) ;
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; fftw_free( in ) ; 
  fftw_cleanup( ) ; 

  // free correlators and momentum list
  free_momcorrs( disp , NSNS , NSNS , NMOM[0] ) ;

  free( NMOM ) ;
  free( (void*)list ) ;

  // free our GAMMAS
  free( GAMMAS ) ;

  // free our spinors
  free( S1 ) ;  free( S1f ) ;
  free( S2 ) ;  free( S2f ) ;

  return FAILURE ;
#endif
}
