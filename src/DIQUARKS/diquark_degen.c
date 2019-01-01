/**
   @file diquark_degen.c
   @brief diquark-diquark contraction for degenerate props
 */
#include "common.h"

#include "basis_conversions.h"    // nrel_rotate_slice()
#include "correlators.h"          // allocate_corrs() && free_corrs()
#include "diquark_contraction.h"  // diquark()
#include "gammas.h"               // make_gammas() && gamma_mmul*
#include "io.h"                   // for read_prop()
#include "progress_bar.h"         // progress_bar()
#include "quark_smear.h"          // sink_smear()
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
  const size_t stride1 = M_CHANNELS ;
  const size_t stride2 = M_CHANNELS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag for if the code messes up
  int error_code = SUCCESS ;

  // loop counters
  size_t i ;

  // gamma LUT
  struct gamma *Cgmu = malloc( stride1 * sizeof( struct gamma ) ) ;
  struct gamma *Cgnu = malloc( stride1 * sizeof( struct gamma ) ) ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  const int sign[ Nprops ] = { +2 } ;
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

  // read in the files
#pragma omp parallel
  {
    // loop counters
    size_t t = 0 , site ;
    
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;

    // smear it if we wish
    sink_smear( M.S , M.S1 , t , CUTINFO , Nprops ) ;

    {
       #pragma omp barrier
    }
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // compute wall sum
      #pragma omp single nowait
      {
	sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
      }

      // assumes all sources are at the same origin, checked in wrap_tetras
      const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ; 

      // read on the master and one slave
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {

	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;

	// loop gamma source
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 * stride2 ; GSGK++ ) {
	  
	  // separate the gamma combinations
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  // perform contraction, result in result
	  M.in[ GSNK + stride2*GSRC ][ site ] = 
	    diquark( SUM_r2[0] , SUM_r2[0] , 
		     Cgmu[ GSRC ] , Cgnu[ GSNK ] , M.GAMMAS[ GAMMA_5 ] ) ;
	}
      }
      // wall-wall contractions
      size_t GSGK ;
      #pragma omp for private(GSGK)
      for( GSGK = 0 ; GSGK < stride1 * stride2 ; GSGK++ ) {
	// separate the gamma combinations
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	M.wwcorr[ GSRC ][ GSNK ].mom[0].C[ tshifted ] = 
	  diquark( M.SUM[0] , M.SUM[1] , Cgmu[ GSRC ] , Cgnu[ GSNK ] ,
		   M.GAMMAS[ GAMMA_5 ] ) ;
      }
      // end of walls

      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ) ;

      // smear the forward prop
      if( t < (LT-1) ) {
	sink_smear( M.Sf , M.S1 , t+1 , CUTINFO , Nprops ) ;
      }

      #pragma omp single
      {
	// copy Sf into S
	copy_props( &M , Nprops ) ;
	
	// status of the computation
	progress_bar( t , LT ) ;
      }
    }
  }

  // write out the diquarks
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;
  
  // memfree sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free C-gammas
  free( Cgmu ) ; free( Cgnu ) ;

  return error_code ;
}

// clean up the number of propagators
#undef Nprops
