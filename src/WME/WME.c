/**
   @file WME.c
   @brief computation of hadronic matrix element using wall sources
 */

#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "contractions.h"      // for the 4 prop contraction
#include "cut_routines.h"      // zero_veclist()
#include "gammas.h"            // gamma matrices
#include "GLU_timer.h"         // print_time() function
#include "io.h"                // read prop
#include "spinor_ops.h"        // spinor multiply
#include "setup.h"             // init_measurements()

// number of propagators
#define Nprops (4)

// ok, brute force this calculation
static double complex
four_quark_trace( const struct spinor SWALL_0 ,
		  const struct spinor DWALL_0 ,
		  const struct spinor SWALL_L_2 ,
		  const struct spinor DWALL_L_2 ,
		  const struct gamma GSRC , 
		  const struct gamma GSNK ,
		  const struct gamma PROJ ,
		  const struct gamma G5 )
{
  // precompute the adjoints, always the "down" quark
  struct spinor anti_DWALL_0 , anti_DWALL_L_2 ;
  full_adj( &anti_DWALL_0 , DWALL_0 , G5 ) ;
  full_adj( &anti_DWALL_L_2 , DWALL_L_2 , G5 ) ;

  // compute
  // trace( SWALL_0 * G5 * anti_DWALL_0 * GSRC * 
  //        SWALL_L_2 * G5 * anti_DWALL_L_2 * GSNK );

  // left and right multiply anti_DWALL_0 by PROJ and GSRC
  gamma_mul_lr( &anti_DWALL_0 , PROJ , GSRC ) ;

  // multiply on the left by SWALL_0
  spinmul_atomic_left( &anti_DWALL_0 , SWALL_0 ) ;

  // left and right multiply anti_DWALL_L_2 by PROJ and GSNK
  gamma_mul_lr( &anti_DWALL_L_2 , PROJ , GSNK ) ;

  // multiply on the left by SWALL_0
  spinmul_atomic_left( &anti_DWALL_L_2 , SWALL_L_2 ) ;

  // and compute the trace of the product of these smaller products
  return bilinear_trace( anti_DWALL_0 , anti_DWALL_L_2 ) ;
}

// attempt to follow UKhadron's implementation where I can.
// Requires one wall at 0 and one well separated at LT/2
// I use the case for BK as an example, e.g. d and s although
// what we put in this is quite general, usually light quark is 
// backward propagating for some reason
int
WME( struct propagator s0 ,
     struct propagator d0 ,
     struct propagator s1 ,
     struct propagator d1 ,
     const struct cut_info CUTINFO ,
     const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices factor of two is because we keep
  // two "terms" of the baryon contraction in "in" and "out"
  const size_t flat_dirac = 2 * stride1 * stride2 ;

  // loop counters
  size_t t ;

  // error code
  int error_code = SUCCESS ;

  // project onto a state : GAMMAS[ 9 ] for projection onto A_t state
  struct gamma PROJ ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { s0 , d0 , s1 , d1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[WME] Failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // pseudoscalar projection state
  PROJ = M.GAMMAS[ GAMMA_5 ] ;

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // read the prop files
#pragma omp parallel
    {
      read_ahead( prop , M.S , &error_code , Nprops ) ;
    }
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // rotate if we must
    rotate_offdiag( M.S , prop , Nprops ) ;

    size_t GSGK ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSGK)
    for( GSGK = 0 ; GSGK < ( NSNS * NSNS ); GSGK++ ) {

      const size_t GSRC = GSGK / ( NSNS ) ;
      const size_t GSNK = GSGK % ( NSNS ) ;
      register double complex trtr = 0.0 ;
      register double complex tr = 0.0 ;
      
      size_t site ;
      for( site = 0 ; site < LCU ; site++ ) {
	// trace-trace component is simple this is projected onto external "PROJ" state
	trtr += ( meson_contract( PROJ , 
				  M.S[1][ site ] , M.GAMMAS[ GSRC ] , 
				  M.S[0][ site ] , M.GAMMAS[ GAMMA_5 ] ) *
		  meson_contract( PROJ , 
				  M.S[3][ site ] , M.GAMMAS[ GSNK ] , 
				  M.S[4][ site ] , M.GAMMAS[ GAMMA_5 ] ) ) ;
	// four quark trace is unpleasant
	tr += four_quark_trace( M.S[0][ site ] , M.S[1][ site ] ,
				M.S[2][ site ] , M.S[3][ site ] ,
				M.GAMMAS[ GSRC ] , M.GAMMAS[ GSNK ] ,
				PROJ , M.GAMMAS[ GAMMA_5 ] ) ;
      }
      // there is probably a factor in this is it 1/2?
      M.corr[ GSRC ][ GSNK ].mom[0].C[ t ] = tr - trtr ;
    }
    
    // tell us how far along we are
    fprintf( stdout , "\r[WME] done %.f %%" , 100.0*(t+1)/(double)(LT) ) ;
    fflush( stdout ) ;
  }
  fprintf( stdout , "\n" ) ;

  // and write out a file
  write_momcorr( outfile , (const struct mcorr**)M.corr ,
		 M.list , stride1 , stride2 , M.nmom , "" ) ;

  // memory deallocation
 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up the number of props
#undef Nprops
