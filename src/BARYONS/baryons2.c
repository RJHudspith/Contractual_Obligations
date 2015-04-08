/**
   @file baryons.c
   @brief baryon contraction code


   This needs some splainin
*/

#include "common.h"

#include "bar_contractions.h"  // cross_color trace and baryon contract
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // gamma_mul_lr()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas() && gamma_mmul()
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()
#include "spinor_ops.h"        // sumprop()

// contracts S2 ( S1 C S1 C ) , odd one out is outside
int
baryons_2diagonal( struct propagator prop1 ,
		   struct propagator prop2 ,
		   const char *outfile )
{
  // gamma matrix storage
  struct gamma *GAMMAS = NULL ;

  // our spinors
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // and the correlators
  struct correlator **Buud_corr = NULL , **Buuu_corr = NULL , **Buds_corr = NULL ;
  struct correlator **Buud_corrWW = NULL , **Buuu_corrWW = NULL , **Buds_corrWW = NULL ;

  // spinor allocations
  if( posix_memalign( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( posix_memalign( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( posix_memalign( (void**)&S2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( posix_memalign( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( prop1.basis == NREL || prop2.basis == NREL ) {
    if( make_gammas( GAMMAS , NREL ) == FAILURE ) {
      goto FREE_FAIL ;
    }
  } else {
   if( make_gammas( GAMMAS , CHIRAL ) == FAILURE ) {
      goto FREE_FAIL ;
    }
  }

  // Define our output correlators, with B_CHANNELS channels and NSNS components
  Buds_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  Buud_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  Buuu_corr = allocate_corrs( B_CHANNELS , NSNS ) ;

  // allocate the walls
  if( prop1.source == WALL ) {
    Buds_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
    Buud_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
    Buuu_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // rotato
    if( prop1.basis == CHIRAL && prop2.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    }
    if( prop2.basis == CHIRAL && prop1.basis == NREL ) {
      nrel_rotate_slice( S2 ) ;
    }				   

    // accumulate wall sum expects both to be walls
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }

    // really need to think about the loop ordering here - J
    int GSRC = 0 ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( L0 - 1 ) ) {
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      #pragma omp single nowait
      {
	if( t < ( L0 - 1 ) ) {
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
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
	const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;

	// accumulate the sums with open dirac indices
	double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
	term[0] = calloc( NSNS , sizeof( double complex ) ) ;
	term[1] = calloc( NSNS , sizeof( double complex ) ) ;
	
	// Wall-Local
	int site ;
	for( site = 0 ; site < LCU ; site++ ) {
	  baryon_contract_site( term , 
				S1[ site ] , S1[ site ] , S2[ site ] ,
				Cgmu , CgmuT ) ;
	}

	// Fill baryon correlator array
	int i ;
	for( i = 0 ; i < NSNS ; i++ ) {
	  Buds_corr[ GSRC ][ i ].C[ t ] = term[0][i] ;
	  Buud_corr[ GSRC ][ i ].C[ t ] = term[0][i] + term[1][i] ;
	  Buuu_corr[ GSRC ][ i ].C[ t ] = 2 * term[0][i] + 4 * term[1][i] ;
	  term[0][ i ] = term[1][ i ] = 0.0 ; // set to zero
	}

	// contract the wall if we desire
	if( prop1.source == WALL ) {
	  baryon_contract_site( term , 
				SUM1 , SUM1 , SUM2 ,
				Cgmu , CgmuT ) ;
	  for( i = 0 ; i < NSNS ; i++ ) {
	    Buds_corrWW[ GSRC ][ i ].C[ t ] = term[0][i] ;
	    Buud_corrWW[ GSRC ][ i ].C[ t ] = term[0][i] + term[1][i] ;
	    Buuu_corrWW[ GSRC ][ i ].C[ t ] = 2 * term[0][i] + 4 * term[1][i] ;
	  }
	}
	// and that is it
	free( term[0] ) ; free( term[1] ) ; free( term ) ;
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
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_baryons( "Baryon: uds-type" , (const struct correlator**)Buds_corr ) ;
  debug_baryons( "Baryon: uud-type" , (const struct correlator**)Buud_corr ) ;
  debug_baryons( "Baryon: uuu-type" , (const struct correlator**)Buuu_corr ) ;
#endif

  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s.uds" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buds_corr , B_CHANNELS , NSNS ) ;
  sprintf( outstr , "%s.uud" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buud_corr , B_CHANNELS , NSNS ) ;
  sprintf( outstr , "%s.uuu" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)Buuu_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;

  // IO for the wall
  if( prop1.source == WALL ) {
    sprintf( outstr , "%s.uds.WW" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)Buds_corrWW , B_CHANNELS , NSNS ) ;
    sprintf( outstr , "%s.uud.WW" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)Buud_corrWW , B_CHANNELS , NSNS ) ;
    sprintf( outstr , "%s.uuu.WW" , outfile ) ;
    write_correlators( outstr , (const struct correlator**)Buuu_corrWW , B_CHANNELS , NSNS ) ;
    free_corrs( Buds_corrWW , B_CHANNELS , NSNS ) ;
    free_corrs( Buud_corrWW , B_CHANNELS , NSNS ) ;
    free_corrs( Buuu_corrWW , B_CHANNELS , NSNS ) ;
  }

  // free stuff
  free( S1 ) ;
  free( S1f ) ;
  free( S2 ) ;
  free( S2f ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

  // failure sink
 FREE_FAIL :

  // free our correlators
  free_corrs( Buds_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buds_corrWW , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corrWW , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corrWW , B_CHANNELS , NSNS ) ;

  // free spinors
  free( S1f ) ;
  free( S1 ) ;
  free( S2f ) ;
  free( S2 ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}
