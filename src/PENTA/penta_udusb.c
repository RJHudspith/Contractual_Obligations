/**
   @file penta_udusb.c
   @brief first pentaquark contraction code

   more or less exactly the same as the tetraquark one
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
#include "penta_contractions.h" // 

// number of propagators
#define Nprops (3)

int
pentaquark_udusb( struct propagator prop1 , // L
		  struct propagator prop2 , // S
		  struct propagator prop3 , // H
		  struct cut_info CUTINFO ,
		  const char *outfile )
{
  // counters
  const size_t stride1 = 2 ;
  const size_t stride2 = PENTA_NOPS*PENTA_NBLOCK*PENTA_NBLOCK ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  uint8_t **loc = malloc( 6561 * sizeof( uint8_t* ) ) ;
  size_t j , k ;
  for( j = 0 ; j < 6561 ; j++ ) {
    loc[ j ] = malloc( 8 * sizeof( uint8_t ) ) ;
    size_t sub = NC , div = 1 ; 
    for( k = 8 ; k > 0 ; k-- ) {
      loc[ j ][ 8 - k ] = (uint8_t)( ( j % sub ) / div ) ;
      sub *= NC ;
      div *= NC ;
    }
  }
  
  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 , prop3 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[TETRA] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // init the parallel region
#pragma omp parallel
  {
    // each thread gets a copy
    double complex **F = malloc( NSNS * sizeof( double complex* ) ) ;
    size_t i ;
    for( i = 0 ; i < NSNS ; i++ ) {
      F[i] = malloc( 6561*sizeof( double complex ) ) ;
    }
    
    // read in the first timeslice
    read_ahead( prop , M.S , &error_code , Nprops ) ;

    #pragma omp barrier
    
    // loop counters
    size_t t ;
  
    // Time slice loop 
    for( t = 0 ; t < LT ; t++ ) {
      
      // if we are doing nonrel-chiral hadrons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ;

      // compute wall sum
      struct spinor SUMbwdH ;
      if( M.is_wall == GLU_TRUE ) {
	#pragma omp single
	{
	  sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
	}
	full_adj( &SUMbwdH , M.SUM[2] , M.GAMMAS[ GAMMA_5 ] ) ;
      }
    
      // assumes all sources are at the same origin, checked in wrap_tetras
      const size_t tshifted = ( t - prop[0].origin[ND-1] + LT ) % LT ;
      
      size_t site ;
      // read on the master and slaves
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      // Loop over spatial volume threads better
#pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	// precompute backward bottom propagator
	struct spinor bwdH ;
	full_adj( &bwdH , M.S[2][ site ] , M.GAMMAS[ GAMMA_5 ] ) ;
	
	// pentaquark contractions stored in result
	double complex *result = malloc( 2 * stride2 * sizeof( double complex ) ) ;
	size_t op ;
	for( op = 0 ; op < 2 * stride2 ; op++ ) {
	  result[ op ] = 0.0 ;
	}

	const struct spinor SUM0_r2 = sum_spatial_sep( M , site , 0 ) ;
	const struct spinor SUM1_r2 = sum_spatial_sep( M , site , 1 ) ;

	// perform contraction, result in result
	pentas( result , F , SUM0_r2 , SUM1_r2 , bwdH , M.GAMMAS ,
		(const uint8_t**)loc ) ;
	
	// put contractions into flattend array for FFT
	for( op = 0 ; op < stride2 ; op++ ) {
	  M.in[ op ][ site ] = result[ op ] ;
	  M.in[ op + stride2 ][ site ] = result[ op + stride2 ] ;
	}

	free( result ) ;
      }

      // wall-wall contractions
      if( M.is_wall == GLU_TRUE ) {
	double complex *result =
	  malloc( 2 * stride2 * sizeof( double complex ) )  ;
	size_t k ;
	for( k = 0 ; k < 2 * stride2 ; k++ ) {
	  result[ k ] = 0.0 ;
	}
	// perform contraction, result in result
	pentas( result , F , M.SUM[0] , M.SUM[1] , SUMbwdH , M.GAMMAS ,
		(const uint8_t**)loc ) ;
	// put contractions into final correlator object
	size_t op ;
	for( op = 0 ; op < stride2 ; op++ ) {
	  M.wwcorr[ 0 ][ op ].mom[ 0 ].C[ tshifted ] = result[ op ] ;
	  M.wwcorr[ 1 ][ op ].mom[ 0 ].C[ tshifted ] = result[ op + stride2 ] ;
	}
	free( result ) ;
      }
      // end of wall-wall stuff
      
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

    for( i = 0 ; i < NSNS ; i++ ) {
      free( F[ i ] ) ;
    }
    free( F ) ;
  }

  // skip writing the files out if we fucked up
  if( error_code == FAILURE ) goto memfree ;

  // write out the penta wall-local and maybe wall-wall
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

// clean up number of props
#undef Nprops