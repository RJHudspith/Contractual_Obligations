/**
   @file mesons.c
   @brief dispersion relation computation for mesons
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gt_Gdag_gt()
#include "io.h"                // read_prop
#include "progress_bar.h"      // progress_bar()
#include "quark_smear.h"       // sink_smear()
#include "setup.h"             // free_ffts() ..
#include "spinor_ops.h"        // sumprop()

// number of propagators
#define Nprops (1)

// computes flavour-diagonal correlators
int
mesons_diagonal( struct propagator prop1 ,
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

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  // phases cancel for mesons of same propagator
  const int sign[ Nprops ] = { 0 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[MESONS] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }  

  // initialise the parallel region
#pragma omp parallel
  {
    // loop counters
    size_t t = 0 , site ;
    
    // initially read in a timeslice
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;

    // smear it if we wish
    sink_smear( M.S , M.S1 , t , CUTINFO , Nprops ) ;

    {
       #pragma omp barrier
    }
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // compute wall-wall sum
      #pragma omp single nowait
      {
	sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
      }
    
      // for multiple time sources
      const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ;
      
      // master-slave the IO and perform each FFT (if available) in parallel
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {

	// sum over possible spatial extensions on the sink side
	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;
	
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  #ifdef TWOPOINT_FILTER
	  if( !filter[ GSRC ][ GSNK ] ) continue ;
	  #endif
	  
	  const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
							 M.GAMMAS[ GAMMA_T ] ) ;
	  // loop spatial hypercube
	  // contract with the summed spinor
	  M.in[ GSGK ][ site ] = 
	    meson_contract( gt_GSNKdag_gt    , SUM_r2[0]  , 
			    M.GAMMAS[ GSRC ] , SUM_r2[0] ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
	// correlator computed just out of the summed walls
      }  
      // end of loop on sites
      size_t GSGK ;
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	#ifdef TWOPOINT_FILTER
	if( !filter[ GSRC ][ GSNK ] ) continue ;
        #endif
	
	const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
						       M.GAMMAS[ GAMMA_T ] ) ;
	M.wwcorr[ GSRC ][ GSNK ].mom[0].C[ tshifted ] =	\
	  meson_contract( gt_GSNKdag_gt  , M.SUM[0] , 
			  M.GAMMAS[ GSRC ] , M.SUM[0] ,
			  M.GAMMAS[ GAMMA_5 ] ) ;
      }

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

  if( error_code == FAILURE ) goto memfree ;
  
  // write out the ND-1 momentum-injected correlator and the wall
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;
  
  return error_code ;
}

// clean up number of props
#undef Nprops
