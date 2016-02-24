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
#include "setup.h"                // *_measurements()
#include "spinor_ops.h"           // sumprop()

// number of propagators for this code
#define Nprops (1) 

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
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag for if the code messes up
  int error_code = SUCCESS ;

  // loop counters
  size_t i , t , site ;

  // gamma LUT
  struct gamma *Cgmu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;
  struct gamma *Cgnu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;

  // wall sums
  struct spinor SUM1 ;

  // initialise our measurement struct
  const struct propagator prop[ Nprops ] = { prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gammas we use
  for( i = 0 ; i < B_CHANNELS ; i++ ) {
    Cgmu[ i ] = CGmu( M.GAMMAS[ i ] , M.GAMMAS ) ;
    Cgnu[ i ] = gt_Gdag_gt( Cgmu[i] , M.GAMMAS[ GAMMA_T ] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , M.S[0] ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall sum
    if( M.is_wall == GLU_TRUE ) {
      sumprop( &SUM1 , M.S[0] ) ;
    }

    // assumes all sources are at the same origin, checked in wrap_tetras
    const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 

    // strange memory access pattern threads better than what was here before
    #pragma omp parallel
    {
      // read on the master and one slave
      if( t < LT-1 ) {
        #pragma omp master
	{
	  if( read_prop( prop1 , M.Sf[0] ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// loop gamma source
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 * stride2 ; GSGK++ ) {
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  // perform contraction, result in result
	  M.in[ GSNK + stride2*GSRC ][ site ] = 
	    diquark( M.S[0][ site ] , M.S[0][ site ] , 
		     Cgmu[ GSRC ] , Cgnu[ GSNK ] ) ;
	}
      }
      // wall-wall contractions
      if( prop1.source == WALL ) {
	// loop gamma source
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 * stride2 ; GSGK++ ) {
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  M.wwcorr[ GSRC ][ GSNK ].mom[0].C[ tshifted ] = 
	    diquark( SUM1 , SUM1 , Cgmu[ GSRC ] , Cgnu[ GSNK ] ) ;
	}
      }
      // end of walls
    }

    // compute the contracted correlator
    compute_correlator( &M , stride1 , stride2 , tshifted ) ;

    // if we error we leave
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &M.S[0][i] , &M.Sf[0][i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf( "\r[DIQUARK] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)M.corr ,
		 M.list , stride1 , stride2 , M.nmom , "" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr ,
		   M.wwlist , stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

  // failure sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free C-gammas
  free( Cgmu ) ; free( Cgnu ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}

// clean up the number of propagators
#undef Nprops
