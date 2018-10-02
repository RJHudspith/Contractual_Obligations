/**
   @file tetra_udcb.c
   @brief tetraquark contraction code
*/

#include "common.h"

#include "basis_conversions.h"  // nrel_rotate_slice()
#include "contractions.h"       // gamma_mul_lr()
#include "correlators.h"        // allocate_corrs() && free_corrs()
#include "cut_routines.h"       // veclist
#include "gammas.h"             // make_gammas() && gamma_mmul*
#include "io.h"                 // for read_prop()
#include "progress_bar.h"       // progress_bar()
#include "setup.h"              // compute_correlator() ..
#include "spinor_ops.h"         // sumprop()
#include "tetra_contractions.h" // diquark_diquark()

// number of propagators
#define Nprops (3)

// tetraquark candidate L1 L1 \bar{H1} \bar{H2} so three propagators
int
tetraquark_udcb( struct propagator prop1 , // L1
		 struct propagator prop2 , // H1
		 struct propagator prop3 , // H2
		 struct cut_info CUTINFO ,
		 const char *outfile )
{
  // counters
  const size_t stride1 = TETRA_NOPS ;
  const size_t stride2 = ND-1 ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 , prop3 } ;
  const int sign[ Nprops ] = { +2 , -1 , -1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[TETRA] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // init the parallel region
#pragma omp parallel
  {
    // loop counters
    size_t t = 0 ;
    
    // read in the first timeslice
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;

    {
       #pragma omp barrier
    }
    
    // Time slice loop 
    for( t = 0 ; t < LT ; t++ ) {
      
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
      // read on the master and slaves
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// tetraquark contractions stored in result
	double complex result[ stride1 ] ;
	size_t GSRC , op ;
	for( op = 0 ; op < stride1 ; op++ ) {
	  result[ op ] = 0.0 ;
	}

	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;

	// precompute backward bottom propagator
	struct spinor bwdH1_r2 , bwdH2_r2 ;
	full_adj( &bwdH1_r2 , SUM_r2[1] , M.GAMMAS[ GAMMA_5 ] ) ;
	full_adj( &bwdH2_r2 , SUM_r2[2] , M.GAMMAS[ GAMMA_5 ] ) ;
	
	// loop gamma source
	for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	  // perform contraction, result in result
	  tetras( result , SUM_r2[0] , SUM_r2[0] , bwdH1_r2 , bwdH2_r2 ,
		  M.GAMMAS , GSRC , GLU_TRUE , GLU_FALSE ) ;
	  // put contractions into flattend array for FFT
	  for( op = 0 ; op < stride1 ; op++ ) {
	    M.in[ GSRC + op * stride2 ][ site ] = result[ op ] ;
	  }
	}
      }
      // wall-wall contractions
      struct spinor SUMbwdH1 , SUMbwdH2 ;
      full_adj( &SUMbwdH1 , M.SUM[1] , M.GAMMAS[ GAMMA_5 ] ) ;
      full_adj( &SUMbwdH2 , M.SUM[2] , M.GAMMAS[ GAMMA_5 ] ) ;
      size_t GSRC  ;
      #pragma omp for private(GSRC)
      for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	double complex result[ stride1 ] ;
	size_t op ;
	for( op = 0 ; op < stride1 ; op++ ) {
	  result[ op ] = 0.0 ;
	}
	// perform contraction, result in result
	tetras( result , M.SUM[0] , M.SUM[0] , SUMbwdH1 , SUMbwdH2 ,
		M.GAMMAS , GSRC , GLU_TRUE , GLU_FALSE ) ;
	// put contractions into final correlator object
	for( op = 0 ; op < stride1 ; op++ ) {
	  M.wwcorr[ op ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[ op ] ;
	}
	///
      }
      
      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ) ;
      
      #pragma omp single
      {
	// copy over the propagators
	copy_props( &M , Nprops ) ;
	
	// status of the computation
	progress_bar( t , LT ) ;
      }
    }
  }

  // skip writing out files if we fucked up
  if( error_code == FAILURE ) goto memfree ;

  // write out the tetra wall-local and maybe wall-wall
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;

  // memfree sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
