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
#include "read_propheader.h"   // (re)read the header
#include "spinor_ops.h"        // spinor multiply

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
     const char *outfile )
{
  // usual cruft
  struct spinor *SWALL_0 = NULL , *SWALL_L_2 = NULL ;
  struct spinor *DWALL_0 = NULL , *DWALL_L_2 = NULL ;

  // allocate the basis
  struct gamma *GAMMAS = NULL ;

  // momentum list stuff
  int *NMOM = NULL ;
  struct veclist *list = NULL ;

  // data structure for holding the contractions
  struct mcorr **corr = NULL ;

  // flag for whether we switch basis
  GLU_bool NREL_FLAG = GLU_FALSE ;

  // allocate our four spinors expecting them to be at 0 and L/2
  if( corr_malloc( (void**)&SWALL_0 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&DWALL_0 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&SWALL_L_2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }
  if( corr_malloc( (void**)&DWALL_L_2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto free_failure ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( s0.basis == NREL || d0.basis == NREL ||
      s1.basis == NREL || d1.basis == NREL ) {
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      goto free_failure ;
    }
    NREL_FLAG = GLU_TRUE ;
  } else {
    if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      goto free_failure ;
    }
  }

  // pseudoscalar projection state
  const struct gamma PROJ = GAMMAS[ GAMMA_5 ] ; // GAMMAS[ 9 ] for projection onto A_t state

  // create a ( 0 , 0 , 0 ) vector list
  NMOM = malloc( sizeof( int ) ) ;
  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

  // allocate our corrs
  corr = allocate_momcorrs( NSNS , NSNS , NMOM[0] ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( s0 , SWALL_0 ) == FAILURE || 
	read_prop( d0 , DWALL_0 ) == FAILURE ||
	read_prop( s1 , SWALL_L_2 ) == FAILURE ||
	read_prop( d1 , DWALL_L_2 ) == FAILURE ) {
      goto free_failure ;
    }

    // if any of these are a mix of chiral && nrel we rotate all to NREL
    if( NREL_FLAG == GLU_TRUE ) {
      if( s0.basis == CHIRAL ) nrel_rotate_slice( SWALL_0 ) ;
      if( d0.basis == CHIRAL ) nrel_rotate_slice( DWALL_0 ) ;
      if( s1.basis == CHIRAL ) nrel_rotate_slice( SWALL_L_2 ) ;
      if( d1.basis == CHIRAL ) nrel_rotate_slice( DWALL_L_2 ) ;
    }

    int GSGK ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSGK)
    for( GSGK = 0 ; GSGK < ( NSNS * NSNS ); GSGK++ ) {

      const int GSRC = GSGK / ( NSNS ) ;
      const int GSNK = GSGK % ( NSNS ) ;
      register double complex trtr = 0.0 ;
      register double complex tr = 0.0 ;
      
      int site ;
      for( site = 0 ; site < LCU ; site++ ) {
	// trace-trace component is simple this is projected onto external "PROJ" state
	trtr += ( meson_contract( PROJ , DWALL_0[ site ] ,
				  GAMMAS[ GSRC ] , SWALL_0[ site ] , 
				  GAMMAS[ GAMMA_5 ] ) *
		  meson_contract( PROJ , DWALL_L_2[ site ] ,
				  GAMMAS[ GSNK ] , SWALL_L_2[ site ] ,
				  GAMMAS[ GAMMA_5 ] ) ) ;
	// four quark trace is unpleasant
	tr += four_quark_trace( SWALL_0[ site ] , DWALL_0[ site ] ,
				SWALL_L_2[ site ] , DWALL_L_2[ site ] ,
				GAMMAS[ GSRC ] , GAMMAS[ GSNK ] ,
				PROJ , GAMMAS[ GAMMA_5 ] ) ;
      }
      // there is probably a factor in this is it 1/2?
      corr[ GSRC ][ GSNK ].mom[0].C[ t ] = tr - trtr ;
    }
    
    // tell us how far along we are
    printf("\r[WME] done %.f %%",100.0*(t+1)/(double)(L0) ) ;
    fflush( stdout ) ;
  }
  printf("\n") ;
#ifdef DEBUG
  debug_mesons( "[WME]" , (const struct correlator**)corr ) ;
#endif

  // and write out a file
  write_momcorr( outfile , (const struct mcorr**)corr ,
		 list , NSNS , NSNS , NMOM ) ;

  // free our correlator measurement
  free_momcorrs( corr , NSNS , NSNS , NMOM[0] ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free momentum list stuff
  free( NMOM ) ; free( (void*)list ) ;

  // free the space
  free( DWALL_0 ) ; free( DWALL_L_2 ) ;
  free( SWALL_0 ) ; free( SWALL_L_2 ) ;

  // reread headers 
  rewind( s0.file ) ; read_propheader( &s0 ) ;
  rewind( d0.file ) ; read_propheader( &d0 ) ;
  rewind( s1.file ) ; read_propheader( &s1 ) ;
  rewind( d1.file ) ; read_propheader( &d1 ) ;

  // tell us how long it all took, my guess is a long time
  print_time( ) ;

  return SUCCESS ;

 free_failure :

  // free our correlator measurement
  if( NMOM != NULL ) {
    free_momcorrs( corr , NSNS , NSNS , NMOM[0] ) ;
  }

  // free momentum list stuff
  free( NMOM ) ; free( (void*)list ) ;

  // free our gamma matrices
  free( GAMMAS ) ;

  // free the space
  free( DWALL_0 ) ; free( DWALL_L_2 ) ;
  free( SWALL_0 ) ; free( SWALL_L_2 ) ;

  return FAILURE ;
}
     
