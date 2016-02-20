/**
   @file setup.c
   @brief common setup/free code goes here
 */
#include "common.h"

#include "cut_routines.h"         // veclist

// compute the momentum-projected correlation function
int
compute_correlator( struct mcorr **corr , 
		    const double complex **in , 
		    const double complex **out , 
		    const struct veclist *list ,
		    const int *NMOM ,
		    const void *forward , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted )
{
#ifdef HAVE_FFTW3_H
  fftw_plan *fwd = (fftw_plan*)forward ;
#endif
  // momentum projection
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < stride1 ; i++ ) {
    size_t j , p ;
    for( j = 0 ; j < stride2 ; j++ ) {
      #ifdef HAVE_FFTW3_H
      fftw_execute( fwd[ j + stride2*i ] ) ;
      for( p = 0 ; p < NMOM[0] ; p++ ) {
	corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	  out[ j + stride2*i ][ list[ p ].idx ] ;
      }
      #else
      register double complex sum = 0.0 ;
      for( p = 0 ; p < LCU ; p++ ) {
	sum += in[ j + stride2*i ][ p ] ;
      }
      corr[ i ][ j ].mom[ 0 ].C[ tshifted ] = sum ;
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

