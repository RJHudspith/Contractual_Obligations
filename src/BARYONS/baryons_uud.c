/**
   @file baryons_uud.c
   @brief two flavour baryon contraction code
*/
#include "common.h"

#include "bar_contractions.h"  // baryon_contract_site_mom()
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // gamma_mul_lr()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas() && gamma_mmul*
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "plan_ffts.h"         // create_plans_DFT() 
#include "read_propheader.h"   // for read_propheader()
#include "setup.h"             // free_ffts() ...
#include "spinor_ops.h"        // sumprop()

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
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gammas we use
  for( i = 0 ; i < B_CHANNELS ; i++ ) {
    Cgmu[ i ] = CGmu( M.GAMMAS[ i ] , M.GAMMAS ) ;
    Cgnu[ i ] = gt_Gdag_gt( Cgmu[i] , M.GAMMAS[ GAMMA_T ] ) ;
  }

  // read in the first timeslice
  if( read_ahead( prop , M.S , Nprops ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    rotate_offdiag( M.S , prop , Nprops ) ;

    // accumulate wall sum expects both to be walls
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL ) {
      sumprop( &SUM1 , M.S[0] ) ;
      sumprop( &SUM2 , M.S[1] ) ;
    }

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;

    // strange memory access pattern threads better than what was here before
    size_t site ;
    #pragma omp parallel
    {
      if( t < ( LT - 1 ) ) {
	error_code = read_ahead( prop , M.Sf , Nprops ) ;
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
	baryon_contract_walls( M.wwcorr , SUM1 , SUM1 , 
			       SUM2 , M.GAMMAS , 
			       tshifted , UUD_BARYON ) ;
      }
    }

    baryon_momentum_project( &M , stride1 , stride2 , tshifted , UUD_BARYON ) ;

    // if we error we leave
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &M.S[0][i] , &M.Sf[0][i] , sizeof( struct spinor ) ) ;
      memcpy( &M.S[1][i] , &M.Sf[1][i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

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

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}

// undefine the number of props
#undef Nprops

