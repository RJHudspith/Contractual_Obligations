/**
   @file diquark.c
   @brief diquark-diquark contraction
 */
#include "common.h"

#include "basis_conversions.h"    // nrel_rotate_slice()
#include "correlators.h"          // allocate_corrs() && free_corrs()
#include "diquark_contraction.h"  // diquark()
#include "gammas.h"               // make_gammas() && gamma_mmul*
#include "io.h"                   // for read_prop()
#include "progress_bar.h"         // progress_bar()
#include "read_propheader.h"      // for read_propheader()
#include "setup.h"                // compute_correlator() ..
#include "spinor_ops.h"           // sumprop()

// number of propagators for this code
#define Nprops (2)

int
diquark_offdiag( struct propagator prop1 ,
		 struct propagator prop2 ,
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
  struct gamma *Cgmu = malloc( stride1 * sizeof( struct gamma ) ) ;
  struct gamma *Cgnu = malloc( stride1 * sizeof( struct gamma ) ) ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 } ;
  const int sign[ Nprops ] = { +1 , +1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[DIQUARK] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gammas we use
  for( i = 0 ; i < stride1 ; i++ ) {
    Cgmu[ i ] = CGmu( M.GAMMAS[ i ] , M.GAMMAS ) ;
    Cgnu[ i ] = gt_Gdag_gt( Cgmu[i] , M.GAMMAS[ GAMMA_T ] ) ;
  }

  // read in the first timeslice
  // read in the files
#pragma omp parallel
  {
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;
  }
  if( error_code == FAILURE ) {
    goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    rotate_offdiag( M.S , prop , Nprops ) ;

    // compute wall sum
    if( prop1.source == WALL ) {
      sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
    }

    // assumes all sources are at the same origin, checked in wrap_tetras
    const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 

    // strange memory access pattern threads better than what was here before
    #pragma omp parallel
    {
      // read on the master and one slave
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;

	// loop gammas
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 * stride2 ; GSGK++ ) {
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  // perform contraction, result in result
	  M.in[ GSNK + stride2*GSRC ][ site ] = 
	    diquark( SUM_r2[0] , SUM_r2[1] , 
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
	    diquark( M.SUM[0] , M.SUM[1] , Cgmu[ GSRC ] , Cgnu[ GSNK ] ) ;
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
    copy_props( &M , Nprops ) ;
 
    // status of the computation
    progress_bar( t , LT ) ;
  }

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)M.corr , M.list ,
		 M.sum_twist , stride1 , stride2 , M.nmom , "" ) ;
  // if we have walls we use them
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr , M.wwlist ,
		   M.sum_twist , stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

  // failure sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free C-gammas
  free( Cgmu ) ; free( Cgnu ) ;

  return error_code ;
}

// clean up the number of propagators
#undef Nprops
