/**
   @file WME.c
   @brief computation of hadronic matrix element using wall sources
 */

#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "contractions.h"      // for the 4 prop contraction
#include "gammas.h"            // gamma matrices
#include "read_propheader.h"   // (re)read the header

// a little simplification
static void
rotate_to_NREL( struct spinor *S ) 
{
  int site ;
#pragma omp parallel for private(site) 
  for( site = 0 ; site < LCU ; site++ ) {
    chiral_to_nrel( &S[ site ] ) ;
  }
  return ;
}

// ok brute force this calculation
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
     const char *outfile )
{
  // usual cruft
  // allocate the basis
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  GLU_bool NREL_FLAG = GLU_FALSE ;
  if( s0.basis == NREL || d0.basis == NREL ||
      s1.basis == NREL || d1.basis == NREL ) {
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      free( GAMMAS ) ;
      return FAILURE ;
    }
    NREL_FLAG = GLU_TRUE ;
  } else {
    if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      free( GAMMAS ) ;
      return FAILURE ;
    }
  }

  // data structure for holding the contractions
  struct correlator **corr = (struct correlator**)malloc( NSNS * sizeof( struct correlator* ) ) ;

  allocate_corrs( corr ) ;

  // allocate our four spinors expecting them to be at 0 and L/2
  struct spinor *SWALL_0 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *DWALL_0 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *SWALL_L_2 = calloc( VOL3 , sizeof( struct spinor ) ) ;
  struct spinor *DWALL_L_2 = calloc( VOL3 , sizeof( struct spinor ) ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( s0 , SWALL_0 ) == FAILURE || 
	read_prop( d0 , DWALL_0 ) == FAILURE ||
	read_prop( s1 , SWALL_L_2 ) == FAILURE ||
	read_prop( d1 , DWALL_L_2 ) == FAILURE ) {
      free( SWALL_0 ) ; free( SWALL_L_2 ) ; free( DWALL_0 ) ; 
      free( DWALL_L_2 ) ; free( GAMMAS ) ; free_corrs( corr ) ;
      return FAILURE ;
    }

    // if any of these are a mix of chiral && nrel we rotate all to NREL
    if( NREL_FLAG == GLU_TRUE ) {
      if( s0.basis == CHIRAL ) rotate_to_NREL( SWALL_0 ) ;
      if( d0.basis == CHIRAL ) rotate_to_NREL( DWALL_0 ) ;
      if( s1.basis == CHIRAL ) rotate_to_NREL( SWALL_L_2 ) ;
      if( d1.basis == CHIRAL ) rotate_to_NREL( DWALL_L_2 ) ;
    }

    const struct gamma PROJ = GAMMAS[ GAMMA_5 ] ; // GAMMAS[ GAMMA_5 + 1 ] for projection onto A0 state

    int GSRC ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < NSNS ; GSRC++ ) {

      int GSNK ;
      for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {

	register double complex trtr = 0.0 ;
	register double complex tr = 0.0 ;

	int site ;
	for( site = 0 ; site < LCU ; site++ ) {
	  // trace-trace component is simple this is projected onto external "PROJ" state
	  trtr += ( meson_trace( PROJ , DWALL_0[ site ] ,
				 GAMMAS[ GSRC ] , SWALL_0[ site ] , 
				 GAMMAS[ GAMMA_5 ] ) *
		    meson_trace( PROJ , DWALL_L_2[ site ] ,
				 GAMMAS[ GSNK ] , SWALL_L_2[ site ] ,
				 GAMMAS[ GAMMA_5 ] ) ) ;
	  // four quark trace is unpleasant
	  tr += four_quark_trace( SWALL_0[ site ] , DWALL_0[ site ] ,
				  SWALL_L_2[ site ] , DWALL_L_2[ site ] ,
				  GAMMAS[ GSRC ] , GAMMAS[ GSNK ] ,
				  PROJ , GAMMAS[ GAMMA_5 ] ) ;
	}
	corr[ GSRC ][ GSNK ].C[ t ] = tr - trtr ;
      }
    }
    printf("\r[WME] done %.f %%",(t+1)/((L0)/100.));fflush(stdout) ;
    //
  }
  printf("\n") ;
#ifdef DEBUG
  debug_mesons( "[WME]" , (const struct correlator**)corr ) ;
#endif

  // and write out a file
  write_correlators( outfile , (const struct correlator**)corr ) ;

  // free our correlator measurement
  free_corrs( corr ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free the space
  free( DWALL_0 ) ; free( DWALL_L_2 ) ;
  free( SWALL_0 ) ; free( SWALL_L_2 ) ;

  // reread headers 
  rewind( s0.file ) ; read_check_header( &s0 , GLU_FALSE ) ;
  rewind( d0.file ) ; read_check_header( &d0 , GLU_FALSE ) ;
  rewind( s1.file ) ; read_check_header( &s1 , GLU_FALSE ) ;
  rewind( d1.file ) ; read_check_header( &d1 , GLU_FALSE ) ;

  // tell us how long it all took, my guess is a long time
  print_time( ) ;

  return SUCCESS ;
}
     