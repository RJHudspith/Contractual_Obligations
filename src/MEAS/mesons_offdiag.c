/**
   @file mesons_offdiag.c
   @brief flavour off-diagonal meson contractions

   split from mesons for clarity
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // write_momcorr()
#include "gammas.h"            // gt_Gdag_gt()
#include "io.h"                // read_ahead()
#include "progress_bar.h"      // progress_bar()
#include "setup.h"             // init_measurements()
#include "spinor_ops.h"        // sumwalls()

// number of propagators
#define Nprops (2)

// computes flavour-off diagonal correlators with ND-1 momentum projection
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const struct cut_info CUTINFO ,
		    const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1*stride2 ;

  // error flag for if the code messes up
  int error_code = SUCCESS ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[MESONS] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // open the parallel region
#pragma omp parallel
  {
    // initial read of a timeslice
    read_ahead( prop , M.S , &error_code , Nprops ) ;

    #pragma omp barrier
    
    // loop counters
    size_t t , site ;
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // if we are doing nonrel-chiral mesons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ;

      // compute wall-wall sum
      if( M.is_wall == GLU_TRUE ) {
	#pragma omp single
	{
	  sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
	}
      }
    
      // support for multiple time sources
      const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ;
      
      // master-slave the IO and perform each FFT in parallel
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }

      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {
	
	// sum over possible rs
	const struct spinor SUM0_r2 = sum_spatial_sep( M , site , 0 ) ;
	
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	  const size_t GSRC = GSGK / stride1 ;
	  const size_t GSNK = GSGK % stride2 ;
	  const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
							 M.GAMMAS[ GAMMA_T ] ) ;
	  M.in[ GSGK ][ site ] = 
	    meson_contract( gt_GSNKdag_gt  , M.S[1][ site ] , 
			    M.GAMMAS[ GSRC ] , SUM0_r2 , 
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
      }
      size_t GSGK ;
      #pragma omp for private(GSGK)
      for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
						       M.GAMMAS[ GAMMA_T ] ) ;
	// and contract the walls
	if( M.is_wall == GLU_TRUE ) {
	  M.wwcorr[ GSRC ][ GSNK ].mom[ 0 ].C[ tshifted ] =	\
	    meson_contract( gt_GSNKdag_gt  , M.SUM[1] ,
			    M.GAMMAS[ GSRC ] , M.SUM[0] ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
      }

      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ,
			  CUTINFO.configspace ) ;
      
      #pragma omp single
      {
	// copy the forward timeslice into the active timeslice
	copy_props( &M , Nprops ) ;
	
	// status of the computation
	progress_bar( t , LT ) ;
      }
      // tloop
    }
    // parallel regione
  }

  if( error_code == FAILURE ) goto memfree ;
  
  // write out the ND-1 momentum-injected correlator
  write_momcorr( outfile , (const struct mcorr**)M.corr , M.list , 
		 stride1 , stride2 , M.nmom , "" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr , M.wwlist , 
		   stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

  // memory freeing part
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
