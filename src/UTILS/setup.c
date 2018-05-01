/**
   @file setup.c
   @brief common setup/free code goes here
 */
#include "common.h"

#include "correlators.h"       // for allocate_corrs and free_corrs
#include "cut_routines.h"      // veclist
#include "gammas.h"            // gamma matrices
#include "geometry.h"          // compute_spacing()
#include "plan_ffts.h"         // ND-1 FFTS
#include "setup.h"             // alphabetising
#include "spinor_ops.h"        // spinor_zero_site()

// free our ffts
static int
free_ffts( double complex **in , 
	   double complex **out , 
	   void *forward ,
	   void *backward ,
	   const size_t flat_dirac )
{
  size_t i ;
  // free fftw stuff
#ifdef HAVE_FFTW3_H
  // free the "in" allocation
  if( in != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_free( in[ i ] ) ;
    }
  }
  fftw_free( in ) ;
  if( out != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      fftw_free( out[ i ] ) ;
    }
  }
  fftw_free( out ) ; 
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
#else
  // free the "in" allocation
  if( in != NULL ) {
    for( i = 0 ; i < flat_dirac ; i++ ) {
      free( in[ i ] ) ;
    }
  }
  free( in ) ;
#endif
  return SUCCESS ;
}

// allocate and initialise the usual cruft
static int
init_moms( struct measurements *M ,
	   const struct cut_info CUTINFO )
{  
  M -> nmom = malloc( sizeof( int ) ) ;
  M -> wwnmom = malloc( sizeof( int ) ) ;

  if( M -> is_wall == GLU_TRUE ) {
    M -> list = (struct veclist*)zero_veclist( M -> nmom , ND-1 , 
					       CUTINFO.configspace ) ;
    M -> wwlist = (struct veclist*)zero_veclist( M -> wwnmom , ND-1 , 
						 CUTINFO.configspace ) ;
  } else {
#ifdef HAVE_FFTW3_H
    M -> list = (struct veclist*)compute_veclist( M -> nmom , CUTINFO , ND-1 , 
						  CUTINFO.configspace ) ;
#else
    M -> list = (struct veclist*)zero_veclist( M -> nmom , ND-1 , 
					       CUTINFO.configspace ) ;
#endif
  }

  // initialise spatial summation list
  const struct cut_info ORBITS = \
    { .type = CUTINFO.type ,
      .max_mom = ( M -> is_wall == GLU_TRUE ? CUTINFO.max_r2 : 0 ) ,
      .cyl_width = CUTINFO.cyl_width ,
      .configspace = GLU_TRUE } ;
  
  M -> rlist = compute_veclist( &(M->NR) , ORBITS , ND-1 , GLU_TRUE ) ;
  
  return SUCCESS ;
}

// compute the momentum-projected correlation function
int
compute_correlator( struct measurements *M , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted ,
		    const GLU_bool configspace )
{
  if( configspace == GLU_TRUE ) {
    size_t idx ;
    #pragma omp for nowait private(idx) schedule(dynamic)
    for( idx = 0 ; idx < stride1*stride2 ; idx++ ) {
      const size_t i = idx/stride2 ;
      const size_t j = idx%stride2 ;
      size_t p ;
      for( p = 0 ; p < M -> nmom[ 0 ] ; p++ ) {
	M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	  M -> in[ j + stride2*i ][ M -> list[ p ].idx ] ;
      }
    }
  } else {
    #ifdef HAVE_FFTW3_H
    fftw_plan *fwd = (fftw_plan*)M -> forward ;
    #endif
    // momentum projection
    size_t idx ;
    #pragma omp for nowait private(idx) schedule(dynamic)
    for( idx = 0 ; idx < stride1*stride2 ; idx++ ) {
      const size_t i = idx/stride2 ;
      const size_t j = idx%stride2 ;
      size_t p ;
      #ifdef HAVE_FFTW3_H
      fftw_execute( fwd[ j + stride2*i ] ) ;
      for( p = 0 ; p < M -> nmom[ 0 ] ; p++ ) {
	M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	  M -> out[ j + stride2*i ][ M -> list[ p ].idx ] ;
      }
      #else
      register double complex sum = 0.0 ;
      for( p = 0 ; p < LCU ; p++ ) {
	sum += M -> in[ j + stride2*i ][ p ] ;
      }
      M -> corr[ i ][ j ].mom[ 0 ].C[ tshifted ] = sum ;
      #endif
    }
  }
  return SUCCESS ;
}

// do a time-slice wide copy of our propagators
void
copy_props( struct measurements *M , 
	    const size_t Nprops )
{
  struct spinor *ptr ;
  size_t mu ;
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    ptr = M -> S[ mu ] ;
    M -> S[ mu ] = M -> Sf[ mu ] ;
    M -> Sf[ mu ] = ptr ;
  }
  return ;
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

  // free the summation list
  if( M->rlist != NULL ) {
    free( (void*)M->rlist ) ;
  }

  // free our spinors
  size_t mu ;
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    free( M->S[ mu ] ) ;  
    free( M->Sf[ mu ] ) ;
  }
  free( M->S ) ; 
  free( M->Sf ) ;

  // free the wall sums if they were allocated
  if( M -> SUM != NULL ) {
    free( M -> SUM ) ;
  }

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
  M -> rlist = NULL ;
  M -> GAMMAS = NULL ;
  M -> in = NULL ; M -> out = NULL ;
  M -> forward = NULL ; M -> backward = NULL ;
  M -> S = NULL ; M -> Sf = NULL ;
  M -> SUM = NULL ;

  // allocate S and Sf the forwards prop
  M -> S  = malloc( Nprops * sizeof( struct spinor* ) ) ;
  M -> Sf = malloc( Nprops * sizeof( struct spinor* ) ) ;

  // allocate spinors
  size_t i ;
  for( i = 0 ; i < Nprops ; i++ ) {
    M -> S[i] = M -> Sf[i] = NULL ;
    if( corr_malloc( (void**)&M -> S[ i ]  , ALIGNMENT , LCU * sizeof( struct spinor ) ) != 0 ||
	corr_malloc( (void**)&M -> Sf[ i ] , ALIGNMENT , LCU * sizeof( struct spinor ) ) != 0 ) {
      error_code = FAILURE ; goto end ;
    }
  }

  // are these wall source props
  M -> is_wall = GLU_FALSE ;
  for( i = 0 ; i < Nprops ; i++ ) {
    if( prop[ i ].source == WALL ) {
      M -> is_wall = GLU_TRUE ;
      break ;
    }
  }
  
#ifdef HAVE_FFTW3_H

  // fftw aligns these
  M -> in  = fftw_malloc( flat_dirac * sizeof( double complex* ) ) ;
  M -> out = fftw_malloc( flat_dirac * sizeof( double complex* ) ) ;

  // create the fftw plans
  M -> forward  = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 
  M -> backward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 

  // allocate FFTW storage
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( flat_dirac ) ; i++ ) {
    M -> in[ i ]  = fftw_malloc( LCU * sizeof( double complex ) ) ; 
    M -> out[ i ] = fftw_malloc( LCU * sizeof( double complex ) ) ; 
    // zero the "in" vector just in case
    size_t j ;
    for( j = 0 ; j < LCU ; j++ ) {
      M -> in[ i ][ j ] = 0.0 ; 
    }
  }

  // create spatial volume fftw plans
  create_plans_DFT( M -> forward , M -> backward , 
		    M -> in , M -> out , flat_dirac , ND-1 ) ;
#else

  // allocate fftw stuff
  M -> in  = malloc( flat_dirac * sizeof( double complex* ) ) ;
  for( i = 0 ; i < flat_dirac ; i++ ) {
    M -> in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#endif

  // allocate the wall sums
  if( M -> is_wall == GLU_TRUE ) {
    M -> SUM = malloc( Nprops * sizeof( struct spinor ) ) ;
  }

  // initialise momentum lists
  if( init_moms( M , CUTINFO ) == FAILURE ) {
    error_code = FAILURE ; goto end ;
  }

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

// sum over spatial indices a spinor
struct spinor
sum_spatial_sep( const struct measurements M ,
		 const size_t site1 ,
		 const size_t prop_idx )
{
  size_t r ;
  struct spinor SUM ;
  spinor_zero_site( &SUM ) ;
  for( r = 0 ; r < M.NR ; r++ ) {
    const size_t site2 = compute_spacing( M.rlist[r].MOM , site1 ,
					  ND-1 ) ;
    add_spinors( &SUM , M.S[prop_idx][site2] ) ;
  }
  return SUM ;
}
