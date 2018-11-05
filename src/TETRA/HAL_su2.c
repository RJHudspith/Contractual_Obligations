/**
   @file tetra_degen.c
   @brief tetraquark contraction code for degenerate light quarks
*/
#include "common.h"

#include "corr_malloc.h"
#include "contractions.h"    // full_adj()
#include "correlators.h"     // allocate_corrs() && free_corrs()
#include "HAL_rhorho.h"
#include "io.h"              // for read_prop()
#include "plan_ffts.h"
#include "progress_bar.h"    // progress_bar()
#include "setup.h"           // compute_correlator()
#include "spinor_ops.h"      // sumprop()
#include "rhoeta_contract.h" // rhoeta_contract()

// number of props
#define Nprops (1)

// su2 rhoeta is only the connected diagrams for the
// rho-eta looking fella di-meson:
// ( \bar\psi_a \gamma_5 \psi_a )( \bar\psi_b \gamma_i \psi_b )
int
HAL_su2( struct propagator prop1 ,
	 struct cut_info CUTINFO ,
	 const char *outfile )
{
  // counters
  const size_t stride1 = 1 ;
  const size_t stride2 = 1 ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  const int sign[ Nprops ] = { +4 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[TETRA] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // FFT all that jazz
  //M.in = fftw_malloc( LCU * sizeof( double complex ) ) ;
  M.out[0] = fftw_malloc( LCU * sizeof( double complex ) ) ;

  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , M.in[0] , M.out[0] , ND-1 ) ;

  //printf( "Block alloc?\n" ) ;
  
  // precompute the diquark blocks, called blk. Very big mallocs in here
  struct spinmatrix **blk11 ;
  struct spinmatrix **blk12 ;
  struct spinmatrix **blk21 ;
  struct spinmatrix **blk22 ;
  size_t i ;
  blk11 = malloc( LCU*sizeof(struct spinmatrix*) ) ;
  blk12 = malloc( LCU*sizeof(struct spinmatrix*) ) ;
  blk21 = malloc( LCU*sizeof(struct spinmatrix*) ) ;
  blk22 = malloc( LCU*sizeof(struct spinmatrix*) ) ;
  for( i = 0 ; i < LCU ; i++ ) {
    blk11[i] = malloc( NCNC*NCNC*sizeof(struct spinmatrix) ) ;
    blk12[i] = malloc( NCNC*NCNC*sizeof(struct spinmatrix) ) ;
    blk21[i] = malloc( NCNC*NCNC*sizeof(struct spinmatrix) ) ;
    blk22[i] = malloc( NCNC*NCNC*sizeof(struct spinmatrix) ) ;
  }
  //printf( "Block alloc?\n" ) ;
  
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
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // compute wall sum
      #pragma omp single nowait
      {
	sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
      }

      // assumes all sources are at the same origin, checked in wrap_tetras
      const size_t tshifted = ( t - prop1.origin[ND-1] + LT ) % LT ; 
      
      // strange memory access pattern threads better than what was here before
      size_t site ;
      // read on the master and one slave
      if( t < LT-1 ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }

      HALrhorho_contract( M.in[0] , M.out[0] , forward , backward ,
			  blk11 , blk12 , blk21 , blk22 ,
			  M.S[0] , M.GAMMAS , 0 , 0 ,
			  M.nmom , M.list ) ;
      
      
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

  // skip writing files if we fucked up
  if( error_code == FAILURE ) goto memfree ;
  
  // write out the tetra wall-local and wall-wall
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;

  // memfree sink
 memfree :

  // fftw free some stuff
  fftw_destroy_plan( forward ) ;
  fftw_destroy_plan( backward ) ;

  // free the blocks
  for( i = 0 ; i < LCU ; i++ ) {
    free( blk11[i] ) ;
    free( blk12[i] ) ;
    free( blk21[i] ) ;
    free( blk22[i] ) ;
  }
  free( blk11 ) ;
  free( blk12 ) ;
  free( blk21 ) ;
  free( blk22 ) ;

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up the number of props
#undef Nprops
