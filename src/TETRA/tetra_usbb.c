/**
   @file tetraquark.c
   @brief tetraquark contraction code
*/

#include "common.h"

#include "basis_conversions.h"  // rotate_offdiag()
#include "contractions.h"       // full_adj()
#include "correlators.h"        // allocate_corrs() && free_corrs()
#include "cut_routines.h"       // veclist
#include "gammas.h"             // make_gammas() && gamma_mmul*
#include "io.h"                 // for read_prop()
#include "progress_bar.h"       // progress_bar()
#include "quark_smear.h"        // sink_smear()
#include "setup.h"              // init_measurements() ..
#include "spinor_ops.h"         // sumprop()
#include "tetra_contractions.h" // diquark_diquark()

// number of propagators
#define Nprops (3)

// tetraquark candidate L1 L2 \bar{H} \bar{H} so three propagators
int
tetraquark_usbb( struct propagator prop1 ,
		 struct propagator prop2 ,
		 struct propagator prop3 ,
		 struct cut_info CUTINFO ,
		 const char *outfile )
{
  // counters
  const size_t stride1 = TETRA_NOPS ;
  const size_t stride2 = ND ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 , prop3 } ;
  const int sign[ Nprops ] = { +1 , +1 , -2 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[TETRA] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // read in the first timeslice
#pragma omp parallel
  {
    // loop counters
    size_t t = 0 ;
    
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;

    // smear it if we wish
    sink_smear( M.S , M.S1 , t , CUTINFO , Nprops ) ;
    
    {
       #pragma omp barrier
    }
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {
      
      // if we are doing nonrel-chiral hadrons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ;

      // compute wall sum
      #pragma omp single nowait
      {
	sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
      }

      // assumes all sources are at the same origin, checked in wrap_tetras
      const size_t tshifted = ( t - prop[0].origin[ND-1] + LT ) % LT ; 
      
      // strange memory access pattern threads better than what was here before
      size_t site ;
      // read on the master and one slave
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// tetraquark contractions stored in result
	double complex result[ stride1 ] ;
	size_t op , GSRC ;
	for( op = 0 ; op < stride1 ; op++ ) {
	  result[ op ] = 0.0 ;
	}

	// like the udbb we must spread out all the quarks
	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;
	
	// precompute barred bottom propagator
	struct spinor bwdH_r2 ;
	full_adj( &bwdH_r2 , SUM_r2[2] , M.GAMMAS[ GAMMA_5 ] ) ;
	
	// loop gamma source
	for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	  // perform contraction, result in result
	  tetras( result , SUM_r2[0] , SUM_r2[1] , bwdH_r2 , bwdH_r2 ,
		  M.GAMMAS , GSRC , GLU_FALSE , GLU_TRUE ) ;
	  // put contractions into flattend array for FFT
	  for( op = 0 ; op < stride1 ; op++ ) {
	    M.in[ GSRC + op * stride2 ][ site ] = result[ op ] ;
	  }
	}
      }
      // wall-wall contractions
      struct spinor SUMbwdH ;
      full_adj( &SUMbwdH , M.SUM[2] , M.GAMMAS[ GAMMA_5 ] ) ;
      size_t GSRC  ;
      #pragma omp for private(GSRC)
      for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	double complex result[ stride1 ] ;
	size_t op ;
	for( op = 0 ; op < stride1 ; op++ ) {
	  result[ op ] = 0.0 ;
	}
	// perform contraction, result in result
	tetras( result , M.SUM[0] , M.SUM[1] , SUMbwdH , SUMbwdH ,
		M.GAMMAS , GSRC , GLU_FALSE , GLU_TRUE ) ;
	// put contractions into final correlator object
	for( op = 0 ; op < stride1 ; op++ ) {
	  M.wwcorr[ op ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[ op ] ;
	}
      }
      
      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ) ;

      // smear the forward prop
      if( t < (LT-1) ) {
	sink_smear( M.Sf , M.S1 , t+1 , CUTINFO , Nprops ) ;
      }
      
      #pragma omp single
      {
	// copy over the propagators
	copy_props( &M , Nprops ) ;
	
	// status of the computation
	progress_bar( t , LT ) ;
      }
    }
  }

  if( error_code == FAILURE ) goto memfree ;
  
  // write out the tetra wall-local and maybe wall-wall
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;

  // failure sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
