/**
   @file baryons_uud.c
   @brief two flavour baryon contraction code
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

// number of propagators
#define Nprops (2)

// 2-flavour degenerate baryon contraction
int
baryons_2fdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
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
  struct propagator prop[ Nprops ] = { prop1 , prop2 } ;
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

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    rotate_offdiag( M.S , prop , Nprops ) ;

    // accumulate wall sum expects both to be walls
    if( prop1.source == WALL ) {
      sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
    }

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;

    // strange memory access pattern threads better than what was here before
    size_t site ;
    #pragma omp parallel
    {
      if( t < ( LT - 1 ) ) {
	read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 ; GSGK++ ) {
	  
	  const size_t GSRC = GSGK / B_CHANNELS ;
	  const size_t GSNK = GSGK % B_CHANNELS ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( M.in , M.S[0][ site ] , M.S[0][ site ] , 
				    M.S[1][ site ] , Cgmu[ GSRC ] , 
				    Cgnu[ GSNK ] , GSGK , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop1.source == WALL ) {
	baryon_contract_walls( M.wwcorr , 
			       M.SUM[0] , M.SUM[0] , M.SUM[1] , 
			       Cgmu , Cgnu , tshifted , UUD_BARYON ) ;
      }
    }

    baryon_momentum_project( &M , stride1 , stride2 , tshifted , UUD_BARYON ) ;

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
		 stride1 , stride2 , M.nmom , "uud" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr , M.wwlist , 
		   stride1 , stride2 , M.wwnmom , "uud.ww" ) ;
  }

  // failure sink
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  // free our LUT
  free( Cgmu ) ; free( Cgnu ) ;

  return error_code ;
}

// undefine the number of props
#undef Nprops

