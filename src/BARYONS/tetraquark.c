/**
   @file tetraquark.c
   @brief tetraquark contraction code
*/

#include "common.h"

#include "bar_contractions.h"  // cross_color trace and baryon contract
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // allocate_corrs() && free_corrs()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas() && gamma_mmul
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()

// flavour diagonal tetraquark
int
tetraquark_diagonal( struct propagator prop ,
		     const char *outfile )
{
  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL ;

  // correlators
  struct correlator **Buud_corr = NULL , **Buuu_corr = NULL ;
  struct correlator **Buud_corrWW = NULL , **Buuu_corrWW = NULL ;

  // allocations
  if( posix_memalign( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( posix_memalign( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  // Define our output correlators, with B_CHANNELS channels and NSNS components
  Buud_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  Buuu_corr = allocate_corrs( B_CHANNELS , NSNS ) ;

  // allocate the walls if we are using wall source propagators
  if( prop.source == WALL ) {
    Buud_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
    Buuu_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
  }

  // read in the first timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // accumulate wall sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    // really need to think about the loop ordering here - J
    int GSRC = 0 ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( L0 - 1 ) ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      // parallelise the furthest out loop -> flatten with lvolume?
      // loop B_CHANNELS
      #pragma omp for private(GSRC) schedule(dynamic)
      for( GSRC = 0 ; GSRC < ( B_CHANNELS ) ; GSRC++ ) {
	
	// precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
	const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	// precompute \gamma_t ( Cg_\mu )^{*} \gamma_t -> \Gamma^{T} in note
	const struct gamma CgmuT = CGmuT( GAMMAS[ GSRC ] , GAMMAS ) ;

	// accumulate the sums with open dirac indices
	double complex term1[ NSNS ] = {} ;
	double complex term2[ NSNS ] = {} ;
	
	// Wall-Local
	int site ;
	for( site = 0 ; site < LCU ; site++ ) {
	  baryon_contract_site( term1 , term2 , 
				S1[ site ] , S1[ site ] , S1[ site ] ,
				Cgmu , CgmuT ) ;
	}

	// Fill baryon correlator array
	int i ;
	for( i = 0 ; i < NSNS ; i++ ) {
	  Buud_corr[ GSRC ][ i ].C[ t ] = term1[ i ] + term2[i] ;
	  Buuu_corr[ GSRC ][ i ].C[ t ] = 2 * term1[ i ] + 4 * term2[ i ] ;
	  term1[ i ] = term2[ i ] = 0.0 ; // set to zero
	}

	// contract the wall if we desire
	if( prop.source == WALL ) {
	  baryon_contract_site( term1 , term2 , 
				SUM1 , SUM1 , SUM1 ,
				Cgmu , CgmuT ) ;
	  for( i = 0 ; i < NSNS ; i++ ) {
	    Buud_corrWW[ GSRC ][ i ].C[ t ] = term1[ i ] + term2[i] ;
	    Buuu_corrWW[ GSRC ][ i ].C[ t ] = 2 * term1[ i ] + 4 * term2[ i ] ;
	  }
	}
	// and that is it
      }
    }

    // try this
    if( error_flag == FAILURE ) {
      goto FREE_FAIL ;
    }

    // copy over the propagators
    int i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct site ) ) ;
    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_baryons( "Baryon: uud-type" , (const struct correlator**)Buud_corr ) ;
  debug_baryons( "Baryon: uuu-type" , (const struct correlator**)Buuu_corr ) ;
#endif

  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s.uud" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buud_corr , B_CHANNELS , NSNS ) ;
  sprintf( outstr , "%s.uuu" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buuu_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;

  if( prop.source == WALL ) {
    sprintf( outstr , "%s.uud.WW" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)Buud_corrWW , B_CHANNELS , NSNS ) ;
    sprintf( outstr , "%s.uuu.WW" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)Buuu_corrWW , B_CHANNELS , NSNS ) ;
    free_corrs( Buud_corrWW , B_CHANNELS , NSNS ) ;
    free_corrs( Buuu_corrWW , B_CHANNELS , NSNS ) ;
  }

  // free stuff
  free( S1 ) ;
  free( S1f ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

  // failure sink
 FREE_FAIL :

  // free our correlators
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corrWW , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corrWW , B_CHANNELS , NSNS ) ;

  // free spinors
  free( S1f ) ;
  free( S1 ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}
