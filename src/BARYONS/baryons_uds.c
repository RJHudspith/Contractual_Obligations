/**
   @file baryons_uds.c
   @brief three (different) flavour baryon contraction code
*/
#include "common.h"

#include "bar_contractions.h"  // baryon_contract_site_mom()
#include "basis_conversions.h" // rotate_offdiag()
#include "correlators.h"       // write_momcorr()
#include "gammas.h"            // make_gammas() && gamma_mmul*
#include "io.h"                // for read_prop()
#include "progress_bar.h"      // progress_bar()
#include "read_propheader.h"   // for read_propheader()
#include "setup.h"             // free_ffts()
#include "spinor_ops.h"        // sumwalls()

// number of props
#define Nprops (3)

// no flavours degenerate baryon contraction, is the most memory expensive
int
baryons_3fdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    struct propagator prop3 ,
		    struct cut_info CUTINFO ,
		    const char *outfile )
{
  // counters
  const size_t stride1 = B_CHANNELS * B_CHANNELS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices factor of two is because we keep
  // two "terms" of the baryon contraction in "in" and "out"
  const size_t flat_dirac = 2 * stride1 * stride2 ;

  // output file
  char bar_outfile[ strlen( outfile ) + 4 ] ;
  sprintf( bar_outfile , "%s.uds" , outfile ) ;
  
  // error code
  int error_code = SUCCESS ;

  // gamma LUT
  struct gamma *Cgmu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;
  struct gamma *Cgnu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 , prop2 , prop3 } ;
  const int sign[ Nprops ] = { +1 , +1 , +1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac , sign ) == FAILURE ) {
    fprintf( stderr , "[BARYONS] measurement initialisation failed\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gammas we use
  size_t i ;
  for( i = 0 ; i < B_CHANNELS ; i++ ) {
    Cgmu[ i ] = CGmu( M.GAMMAS[ i ] , M.GAMMAS ) ;
    Cgnu[ i ] = gt_Gdag_gt( Cgmu[i] , M.GAMMAS[ GAMMA_T ] ) ;
  }

  // open parallel region
#pragma omp parallel
  {
    // loop counters
    size_t t = 0 ;
    
    // read in the first timeslice
    read_ahead( prop , M.S , &error_code , Nprops , t ) ;

    // barrier to make sure stuff is read in first
    #pragma omp barrier
    
    // Time slice loop 
    for( t = 0 ; t < LT && error_code == SUCCESS ; t++ ) {

      // if we are doing nonrel-chiral hadrons we switch chiral to nrel
      rotate_offdiag( M.S , prop , Nprops ) ;

      // accumulate wall sum expects both to be walls
      #pragma omp single nowait
      {
	sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
      }

      // multiple time source support
      const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;

      // strange memory access pattern threads better than what was here before
      size_t site ;
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops , t+1 ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {

	// summation of r^2 arrays
	struct spinor SUM_r2[ Nprops ] ;
	sum_spatial_sep( SUM_r2 , M , site ) ;
	
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < ( B_CHANNELS * B_CHANNELS ) ; GSGK++ ) {

	  const size_t GSRC = GSGK / B_CHANNELS ;
	  const size_t GSNK = GSGK % B_CHANNELS ;

	  // Wall-Local
	  baryon_contract_site_mom( M.in , 
				    SUM_r2[0] , SUM_r2[1] , SUM_r2[2] , 
				    Cgmu[ GSRC ] , Cgnu[ GSNK ] , GSGK , 
				    site ) ;
	}
      }
      // loop over open indices performing wall contraction
      baryon_contract_walls( M.corr , 
			     M.SUM[0] , M.SUM[1] , M.SUM[2] , 
			     Cgmu , Cgnu , tshifted , UDS_BARYON ) ;

      // momentum projection 
      baryon_momentum_project( &M , stride1 , stride2 ,
			       tshifted , UDS_BARYON ,
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

  if( error_code == FAILURE ) goto memfree ;
  
  // write out the baryons wall-local and wall-wall  
  write_momcorr_WW( M , bar_outfile , stride1 , stride2 ) ;

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free our LUT
  free( Cgmu ) ; free( Cgnu ) ;

  return error_code ;
}

// undefine the number of props
#undef Nprops
