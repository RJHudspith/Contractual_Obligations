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
#include "setup.h"             // free_ffts() ..
#include "spinor_ops.h"        // sumprop()
#include "spinmatrix_ops.h"    // gammaspinmatrix_trace()

// number of propagators
#define Nprops (1)

// computes flavour-diagonal correlators
int
mesons_diagonal( struct propagator prop1 ,
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

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[MESONS] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // initialise the parallel region
#pragma omp parallel
  {
    // initially read in a timeslice
    read_ahead( prop , M.S , &error_code , Nprops ) ;

    #pragma omp barrier
    
    // loop counters
    size_t t , site ;
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {
      
      // compute wall-wall sum
      if( M.is_wall == GLU_TRUE ) {
	#pragma omp single nowait
	{
	  sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
	}
      }
    
      // for multiple time sources
      const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ;
      
      // master-slave the IO and perform each FFT (if available) in parallel
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {

	// sum over possible spatial extensions on the sink side
	const struct spinor SUM0_r2 = sum_spatial_sep( M , site , 0 ) ;

	size_t GSGK ;
	for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
							 M.GAMMAS[ GAMMA_T ] ) ;
	  // loop spatial hypercube
	  // contract with the summed spinor
	  M.in[ GSGK ][ site ] = 
	    meson_contract( gt_GSNKdag_gt  , M.S[0][ site ] , 
			    M.GAMMAS[ GSRC ] , SUM0_r2 ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}	  
	// correlator computed just out of the summed walls
      }  
      // end of loop on sites
      size_t GSGK ;
      #pragma omp for private(GSGK)
      for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
						       M.GAMMAS[ GAMMA_T ] ) ;
	if( M.is_wall == GLU_TRUE ) {
	  M.wwcorr[ GSRC ][ GSNK ].mom[0].C[ tshifted ] =	\
	    meson_contract( gt_GSNKdag_gt  , M.SUM[0] , 
			    M.GAMMAS[ GSRC ] , M.SUM[0] ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
      }

      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ,
			  CUTINFO.configspace ) ;
      
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
  
  // write out the ND-1 momentum-injected correlator and maybe the wall
  write_momcorr( outfile , (const struct mcorr**)M.corr , 
		 M.list , stride1 , stride2 , M.nmom , "" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr ,
		   M.wwlist , stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
