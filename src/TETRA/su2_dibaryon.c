/**
   @file tetra_degen.c
   @brief tetraquark contraction code for degenerate light quarks
*/
#include "common.h"

#include "basis_conversions.h"     // rotate_offdiag()
#include "correlators.h"           // allocate_corrs() && free_corrs()
#include "cut_routines.h"          // veclist
#include "dibaryon_contractions.h" // dibaryon_contract()
#include "io.h"                    // for read_prop()
#include "progress_bar.h"          // progress_bar()
#include "setup.h"                 // compute_correlator() ..
#include "spinor_ops.h"            // sumprop()

// number of props
#define Nprops (1)

// su2 dibaryon is ( \psi_a C\gamma_i \psi_b )( \psi_a C\gamma_i \psi_b )
// with a sum over gamma index "i"
int
su2_dibaryon( struct propagator prop1 ,
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

      // if we are doing nonrel-chiral mesons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ; 

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
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {

	// have to sum all propagators for the blob sinks - J
	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;
	
	// the trick here is that 0,1 == 1,0 so there is just a factor of 2
	register double complex sum = 0.0 ; 
        sum +=   dibaryon_contract( SUM_r2[0] , M.GAMMAS , 0 , 0 ) ;
	sum += 2*dibaryon_contract( SUM_r2[0] , M.GAMMAS , 0 , 1 ) ;
        sum += 2*dibaryon_contract( SUM_r2[0] , M.GAMMAS , 0 , 2 ) ;
	sum +=   dibaryon_contract( SUM_r2[0] , M.GAMMAS , 1 , 1 ) ;
	sum += 2*dibaryon_contract( SUM_r2[0] , M.GAMMAS , 1 , 2 ) ;
	sum +=   dibaryon_contract( SUM_r2[0] , M.GAMMAS , 2 , 2 ) ;
	M.in[0][ site ] = sum ;
      }

      // have to do wall-wall contraction on a single thread
      #pragma omp single
      {
	register double complex sum = 0.0 ; 
        sum +=   dibaryon_contract( M.SUM[0] , M.GAMMAS , 0 , 0 ) ;
	sum += 2*dibaryon_contract( M.SUM[0] , M.GAMMAS , 0 , 1 ) ;
        sum += 2*dibaryon_contract( M.SUM[0] , M.GAMMAS , 0 , 2 ) ;
	sum +=   dibaryon_contract( M.SUM[0] , M.GAMMAS , 1 , 1 ) ;
	sum += 2*dibaryon_contract( M.SUM[0] , M.GAMMAS , 1 , 2 ) ;
	sum +=   dibaryon_contract( M.SUM[0] , M.GAMMAS , 2 , 2 ) ;
	M.wwcorr[0][0].mom[0].C[ tshifted ] = sum ;
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

  // skip writing files if we fucked up
  if( error_code == FAILURE ) goto memfree ;
  
  // write out the tetra wall-local and wall-wall
  write_momcorr_WW( M , outfile , stride1 , stride2 ) ;

  // memfree sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up the number of props
#undef Nprops
