/**
   @file baryons.c
   @brief baryon contraction code

   TODO :: flatten some part of the src/snk indices and master/slave the IO
   seems tough will need to think about it ...
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
#include "spinor_ops.h"        // sumprop()

// flavour degenerate baryon contraction
int
baryons_diagonal( struct propagator prop ,
		  const char *outfile )
{
  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL ;

  // correlators
  struct correlator **Buud_corr = NULL , **Buuu_corr = NULL , **Buds_corr = NULL ;
  struct correlator **Buud_corrWW = NULL , **Buuu_corrWW = NULL , **Buds_corrWW = NULL ;

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
  Buds_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  Buud_corr = allocate_corrs( B_CHANNELS , NSNS ) ;
  Buuu_corr = allocate_corrs( B_CHANNELS , NSNS ) ;

  // allocate the walls if we are using wall source propagators
  if( prop.source == WALL ) {
    Buds_corrWW = allocate_corrs( B_CHANNELS , NSNS ) ;
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
	const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;

	// accumulate the sums with open dirac indices
	double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
	term[0] = calloc( NSNS , sizeof( double complex ) ) ;
	term[1] = calloc( NSNS , sizeof( double complex ) ) ;

	// Wall-Local
	int site ;
	for( site = 0 ; site < LCU ; site++ ) {
	  baryon_contract_site( term , 
				S1[ site ] , S1[ site ] , S1[ site ] ,
				Cgmu , CgmuT ) ;
	}

	// Fill baryon correlator array
	int i ;
	for( i = 0 ; i < NSNS ; i++ ) {
	  Buds_corr[ GSRC ][ i ].C[ t ] = term[0][i] ;
	  Buud_corr[ GSRC ][ i ].C[ t ] = term[0][i] + term[1][i] ;
	  Buuu_corr[ GSRC ][ i ].C[ t ] = 2 * term[0][i] + 4 * term[1][i] ;
	  term[0][i] = term[1][i] = 0.0 ; // set to zero
	}

	// contract the wall if we desire
	if( prop.source == WALL ) {
	  baryon_contract_site( term ,
				SUM1 , SUM1 , SUM1 ,
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
  free_corrs( Buds_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buud_corr , B_CHANNELS , NSNS ) ;
  free_corrs( Buuu_corr , B_CHANNELS , NSNS ) ;

  if( prop.source == WALL ) {
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

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

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

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}

#if 0
/*
  Omega graveyard, so spooky

  struct spinor CgS, CgS51, CgS52 ;	
  struct spinor DiQ, DiQ51, DiQ52 ;	

  term3 = 0.0 ;
  term4 = 0.0 ;
  term5 = 0.0 ;
  term6 = 0.0 ;

  Omega stuff
  // In case of the Omega ( s ( s Cg5 s )) we need the source and sink separately
  CgS51 = Cgamma_snk( S1[ site ] , GSRC ); -> gamma_mul_l( S1[ site ] , CgmuT ) 
  CgS52 = Cgamma_src( S1[ site ] , GSRC ); -> gamma_mul_r( S1[ site ] , Cgmu ) 

  // There are four more terms for the Omega					
  term3 += baryon_contract( DiQ51, S1[ site ], dirac , dirac , OD1 , OD2 );
  term4 += baryon_contract( DiQ52, CgS52, dirac , OD1 , dirac , OD2 );
  term5 += baryon_contract( DiQ52, CgS52, OD1 , dirac , dirac , OD2 );
  term6 += baryon_contract( DiQ51, S1[ site ], OD1 , dirac , dirac , OD2 );

  Buuu += term1 + term2 + term3 + term4 + term5 + term6;

  // In case of the Omega we need two additional diquarks
  DiQ51 = cross_color_trace( CgS52, CgS51 );
  DiQ52 = cross_color_trace( S1[ site ], CgS51 );
 */
#endif