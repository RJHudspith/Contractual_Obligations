/**
   @file baryons_uuu.c
   @brief baryon contraction code
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

// number of props in this code
#define Nprops (1)

// flavour degenerate baryon contraction
int
baryons_diagonal( struct propagator prop1 ,
		  struct cut_info CUTINFO ,
		  const char *outfile )
{
  // counters
  const size_t stride1 = B_CHANNELS * B_CHANNELS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = 2 * stride1 * stride2 ;

  // error flag
  int error_code = SUCCESS ;

  // loop counters
  size_t i , t ;

  // gamma LUT
  struct gamma *Cgmu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;
  struct gamma *Cgnu = malloc( B_CHANNELS * sizeof( struct gamma ) ) ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[BARYONS] measurement initialisation failed\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gammas we use
  for( i = 0 ; i < B_CHANNELS ; i++ ) {
    Cgmu[ i ] = CGmu( M.GAMMAS[ i ] , M.GAMMAS ) ;
    Cgnu[ i ] = gt_Gdag_gt( Cgmu[i] , M.GAMMAS[ GAMMA_T ] ) ;
  }

  // read in the first timeslice
#pragma omp parallel
  {
    read_ahead( prop , M.S , &error_code , Nprops ) ;
  }
  if( error_code == FAILURE ) {
    goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall-wall sum
    if( M.is_wall == GLU_TRUE ) {
      sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
    }

    // shifted times
    const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;
 
    // strange memory access pattern threads better than what was here before
    size_t site ;
    #pragma omp parallel
    {
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site)
      for( site = 0 ; site < LCU ; site++ ) {

	struct spinor SUM0_r2 = sum_spatial_sep( M , site , 0 ) ;
	
	size_t GSGK ; // combined gamma source and sink indices
	for( GSGK = 0 ; GSGK < stride1 ; GSGK++ ) {

	  // source and sink gammas
	  const size_t GSRC = GSGK / B_CHANNELS ;
	  const size_t GSNK = GSGK % B_CHANNELS ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( M.in ,
				    M.S[0][site] , SUM0_r2 , SUM0_r2 ,
				    Cgmu[ GSRC ] , Cgnu[ GSNK ] , GSGK ,
				    site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( M.is_wall == GLU_TRUE ) {
	baryon_contract_walls( M.wwcorr , 
			       M.SUM[0] , M.SUM[0] , M.SUM[0] , 
			       Cgmu , Cgnu , tshifted , UUU_BARYON ) ;
      }
    }

    // momentum projection 
    baryon_momentum_project( &M , stride1 , stride2 , tshifted , UUU_BARYON ,
			     CUTINFO.configspace ) ;

    // if we error we leave
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    copy_props( &M , Nprops ) ;

    // status of the computation
    progress_bar( t , LT ) ;
  }

  // write out the baryons wall-local and maybe wall-wall
  write_momcorr( outfile , (const struct mcorr**)M.corr , M.list , 
		 stride1 , stride2 , M.nmom , "uuu" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr , M.wwlist , 
		   stride1 , stride2 , M.wwnmom , "uuu.ww" ) ;
  }

  // memory frees
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free our LUT
  free( Cgmu ) ; free( Cgnu ) ;

  return error_code ;
}

// undefine the number of props
#undef Nprops
