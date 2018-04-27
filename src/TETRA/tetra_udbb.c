/**
   @file tetra_degen.c
   @brief tetraquark contraction code for degenerate light quarks
*/

#include "common.h"

#include "basis_conversions.h"  // rotate_offdiag()
#include "contractions.h"       // full_adj()
#include "correlators.h"        // allocate_corrs() && free_corrs()
#include "cut_routines.h"       // veclist
#include "gammas.h"             // make_gammas() && gamma_mmul*
#include "io.h"                 // for read_prop()
#include "progress_bar.h"       // progress_bar()
#include "setup.h"              // compute_correlator() ..
#include "spinor_ops.h"         // sumprop()
#include "tetra_contractions.h" // diquark_diquark()

// number of props
#define Nprops (2)

// tetraquark candidate L1 L1 \bar{H} \bar{H} so two propagators
int
tetraquark_udbb( struct propagator prop1 ,
		 struct propagator prop2 ,
		 struct cut_info CUTINFO ,
		 const char *outfile )
{
  // counters
  const size_t stride1 = TETRA_NOPS ;
  const size_t stride2 = (ND-1) ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[TETRA] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // init the parallel region
  #pragma omp parallel
  {
    // read in the first timeslice
    read_ahead( prop , M.S , &error_code , Nprops ) ;

    #pragma omp barrier
    
    // loop counters
    size_t t ;
  
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // if we are doing nonrel-chiral mesons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ; 

      // compute wall sum
      struct spinor SUMbwdH ;
      if( M.is_wall == GLU_TRUE ) {
	#pragma omp single
	{
	  sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
	}
	full_adj( &SUMbwdH , M.SUM[1] , M.GAMMAS[ GAMMA_5 ] ) ;
      }

      // assumes all sources are at the same origin, checked in wrap_tetras
      const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 
      
      // strange memory access pattern threads better than what was here before
      size_t site ;
      // read on the master and one slave
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// precompute backward bottom propagator using 
	// gamma_5 hermiticity
	struct spinor bwdH ;
	full_adj( &bwdH , M.S[1][ site ] , M.GAMMAS[ GAMMA_5 ] ) ;

	// tetraquark contractions stored in "result"
	double complex result[ stride1 ] ;
	size_t GSRC , op ;
	for( op = 0 ; op < stride1 ; op++ ) {
	  result[ op ] = 0.0 ;
	}
	// for now I keep the anti-Bs at x and vary u and d - J
	const struct spinor SUM0_r2 = sum_spatial_sep( M , site , 0 ) ;
	const struct spinor SUM1_r2 = sum_spatial_sep( M , site , 1 ) ;
	struct spinor bwdH_r2 ;
	full_adj( &bwdH_r2 , SUM1_r2 , M.GAMMAS[ GAMMA_5 ] ) ;
	
	// loop gamma source
	for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	  // perform contraction, result in result
	  tetras( result , SUM0_r2 , SUM0_r2 , bwdH_r2 , bwdH ,
		  M.GAMMAS , GSRC , GLU_TRUE , GLU_TRUE ) ;
	  // put contractions into flattend array for FFT
	  for( op = 0 ; op < stride1 ; op++ ) {
	    M.in[ GSRC + op * stride2 ][ site ] = result[ op ] ;
	  }
	}
      }
      // wall-wall contractions
      if( M.is_wall == GLU_TRUE ) {
	size_t GSRC ;
        #pragma omp for private(GSRC)
	for( GSRC = 0 ; GSRC < stride2 ; GSRC++ ) {
	  double complex result[ stride1 ] ;
	  size_t op ;
	  for( op = 0 ; op < stride1 ; op++ ) {
	    result[ op ] = 0.0 ;
	  }
	  // perform contraction, result in result
	  tetras( result , M.SUM[0] , M.SUM[0] , SUMbwdH , SUMbwdH ,
		  M.GAMMAS , GSRC , GLU_TRUE , GLU_TRUE ) ;
	  // put contractions into final correlator object
	  for( op = 0 ; op < stride1 ; op++ ) {
	    M.wwcorr[ op ][ GSRC ].mom[ 0 ].C[ tshifted ] = result[ op ] ;
	  }
	}
      }

      // compute the contracted correlator
      compute_correlator( &M , stride1 , stride2 , tshifted ,
			  CUTINFO.configspace ) ;
      
      #pragma omp single
      {
	// copy over the propagators
	copy_props( &M , Nprops ) ;
	
	// status of the computation
	progress_bar( t , LT ) ;
      }
    }
  }

  // skip writing files if we fucked up
  if( error_code == FAILURE ) goto memfree ;
  
  // write out the tetra wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)M.corr ,
		 M.list , stride1 , stride2 , M.nmom , "" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr ,
		   M.wwlist , stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

  // memfree sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up the number of props
#undef Nprops
