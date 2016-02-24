/**
   @file setup.c
   @brief common setup/free code goes here
 */
#include "common.h"

#include "correlators.h"       // for allocate_corrs and free_corrs
#include "cut_routines.h"      // veclist
#include "gammas.h"            // gamma matrices
#include "plan_ffts.h"         // ND-1 FFTS
#include "setup.h"             // alphabetising

// compute the momentum-projected correlation function
int
compute_correlator( struct measurements *M , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted )
{
#ifdef HAVE_FFTW3_H
  fftw_plan *fwd = (fftw_plan*)M -> forward ;
#endif
  // momentum projection
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < stride1 ; i++ ) {
    size_t j , p ;
    for( j = 0 ; j < stride2 ; j++ ) {
      #ifdef HAVE_FFTW3_H
      fftw_execute( fwd[ j + stride2*i ] ) ;
      for( p = 0 ; p < M -> nmom[ 0 ] ; p++ ) {
	M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	  M -> out[ j + stride2*i ][ M -> list[ p ].idx ] ;
      }
      #else
      register double complex sum = 0.0 ;
      for( p = 0 ; p < LCU ; p++ ) {
	sum += M - > in[ j + stride2*i ][ p ] ;
      }
      M -> corr[ i ][ j ].mom[ 0 ].C[ tshifted ] = sum ;
      #endif
    }
  }
  return SUCCESS ;
}

// free our ffts
int
free_ffts( double complex **in , 
	   double complex **out , 
	   void *forward ,
	   void *backward ,
	   const size_t flat_dirac )
{
  size_t i ;
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
  free( out ) ; 
  if( forward != NULL ) {
    fftw_plan *forwd = (fftw_plan*)forward ;
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_destroy_plan( forwd[i] ) ;
    }
  }
  free( forward )  ; 
  if( backward != NULL ) {
    fftw_plan *bckwd = (fftw_plan*)backward ;
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_destroy_plan( bckwd[i] ) ;
    }
  }
  free( backward ) ; 
  fftw_cleanup( ) ; 
#endif
  return SUCCESS ;
}

// free our measurement struct
void
free_measurements( struct measurements *M ,
		   const size_t Nprops ,
		   const size_t stride1 , 
		   const size_t stride2 , 
		   const size_t flat_dirac )
{
  // free correlators and momentum list
  if( M -> nmom != NULL ) {
    free_momcorrs( M -> corr , stride1 , stride2 , M -> nmom[0] ) ;
    if( M -> is_wall == GLU_TRUE ) {
      free_momcorrs( M -> wwcorr , stride1 , stride2 , M -> wwnmom[0] ) ;
    }
  }

  // free our ffts
  free_ffts( M->in , M->out , M->forward , M->backward , flat_dirac ) ;

  // free momenta lists
  if( M-> nmom != NULL ) {
    free( M->nmom ) ; 
    free( (void*)M->list ) ;
  }
  if( M -> wwnmom != NULL ) {
    free( M->wwnmom ) ; 
    free( (void*)M->wwlist ) ;
  }

  // free our GAMMAS
  free( M->GAMMAS ) ;

  // free our spinors
  size_t mu ;
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    free( M->S[ mu ] ) ;  
    free( M->Sf[ mu ] ) ;
  }
  free( M->S ) ; 
  free( M->Sf ) ;

  return ;
}

// initialise our measurement struct
int
init_measurements( struct measurements *M ,
		   const struct propagator *prop ,
		   const size_t Nprops ,
		   const struct cut_info CUTINFO ,
		   const size_t stride1 ,
		   const size_t stride2 ,
		   const size_t flat_dirac ) 
{
  // error code
  int error_code = SUCCESS ;

  // nullify everything
  M -> nmom = NULL ; M -> wwnmom = NULL ;
  M -> list = NULL ; M -> wwlist = NULL ;
  M -> corr = NULL ; M -> wwcorr = NULL ;
  M -> GAMMAS = NULL ;
  M -> in = NULL ; M -> out = NULL ;
  M -> forward = NULL ; M -> backward = NULL ;
  M -> S = NULL ; M -> Sf = NULL ;

  // allocate S and Sf the forwards prop
  M -> S  = malloc( Nprops * sizeof( struct spinor* ) ) ;
  M -> Sf = malloc( Nprops * sizeof( struct spinor* ) ) ;

  // allocate spinors
  size_t i ;
  for( i = 0 ; i < Nprops ; i++ ) {
    M -> S[i] = M -> Sf[i] = NULL ;
    if( corr_malloc( (void**)&M -> S[ i ]  , 16 , LCU * sizeof( struct spinor ) ) != 0 ||
	corr_malloc( (void**)&M -> Sf[ i ] , 16 , LCU * sizeof( struct spinor ) ) != 0 ) {
      error_code = FAILURE ; goto end ;
    }
  }

  // allocate fftw stuff
  M -> in  = (double complex**)malloc( flat_dirac * sizeof( double complex* ) ) ;
  for( i = 0 ; i < flat_dirac ; i++ ) {
    M -> in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H
  M -> out = (double complex**)malloc( flat_dirac * sizeof( double complex* ) ) ;

  M -> forward  = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 
  M -> backward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( flat_dirac ) ; i++ ) {
    M -> out[ i ] = ( double complex* )malloc( LCU * sizeof( double complex ) ) ; 
  }

  // create spatial volume fftw plans
  create_plans_DFT( M -> forward , M -> backward , 
		    M -> in , M -> out , flat_dirac , ND-1 ) ;
#endif

  // are these wall source props
  M -> is_wall = GLU_FALSE ;
  for( i = 0 ; i < Nprops ; i++ ) {
    if( prop[ i ].source == WALL ) {
      M -> is_wall = GLU_TRUE ;
      break ;
    }
  }

  // initialise momentum lists
  init_moms( &M -> nmom , &M -> wwnmom , 
	     &M -> list , &M -> wwlist , 
	     CUTINFO , M -> is_wall ) ;

  // allocate correlators
  M -> corr = allocate_momcorrs( stride1 , stride2 , M -> nmom[0] ) ;
  if( M -> is_wall == GLU_TRUE ) {
    M -> wwcorr = allocate_momcorrs( stride1 , stride2 , M -> wwnmom[0] ) ;
  }

  // precompute the gamma basis
  M -> GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( setup_gamma( M -> GAMMAS , prop , Nprops ) == FAILURE ) {
    error_code = FAILURE ; goto end ;
  }

 end :
  return error_code ;
}

// allocate and initialise the usual cruft
int
init_moms( int **NMOM , 
	   int **wwNMOM ,
	   struct veclist **list ,
	   struct veclist **wwlist ,
	   const struct cut_info CUTINFO , 
	   const GLU_bool is_wall )
{
  *(NMOM) = malloc( sizeof( int ) ) ;
  *(wwNMOM) = malloc( sizeof( int ) ) ;
#ifdef HAVE_FFTW3_H
  *(list) = (struct veclist*)compute_veclist( *NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;
#else
  *(list) = (struct veclist*)zero_veclist( *NMOM , ND-1 , GLU_FALSE ) ;
#endif
  if( is_wall == GLU_TRUE ) {
    *(wwlist) = (struct veclist*)zero_veclist( *wwNMOM , ND-1 , GLU_FALSE ) ;
  }
  return SUCCESS ;
}
