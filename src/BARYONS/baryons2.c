/**
   @file baryons2.c
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
  const size_t flat_dirac = 2*stride1*stride2 ;

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // momentum stuff
  int *NMOM = NULL , *wwNMOM = NULL ;
  struct veclist *list = NULL , *wwlist = NULL ;

  // fftw temporaries
  double complex **in = NULL , **out = NULL ;

  // correlators
  struct mcorr **Buud_corr = NULL , **Buud_corrWW = NULL ;

#ifdef HAVE_FFTW3_H
  fftw_plan *forward = NULL , *backward = NULL ; 
#else
  int *forward = NULL , *backward = NULL ;
#endif

  // loops
  size_t i , t ;

  // error code
  int error_code = SUCCESS ;

  // allocations
  if( corr_malloc( (void**)&S1  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 || 
      corr_malloc( (void**)&S2  , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    error_code = FAILURE ; goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( setup_gamma_2( GAMMAS , prop1.basis , prop2.basis ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // allocate results matrix
  in = malloc( ( flat_dirac ) * sizeof( double complex* ) ) ;
  for( i = 0 ; i < ( flat_dirac ) ; i++ ) {
    in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H

  out = malloc( ( flat_dirac ) * sizeof( double complex* ) ) ;
  for( i = 0 ; i < ( flat_dirac ) ; i++ ) {
    out[ i ] = fftw_malloc( LCU * sizeof( double complex ) ) ;
  }

  forward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( flat_dirac * sizeof( fftw_plan ) ) ;

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , flat_dirac , ND-1 ) ;

#endif

  // initialise momentum lists
  init_moms( &NMOM , &wwNMOM , &list , &wwlist , CUTINFO , 
	     prop1.source == WALL ? GLU_TRUE : GLU_FALSE ) ;

  // Define our output correlators, with B_CHANNELS channels and NSNS components
  Buud_corr = allocate_momcorrs( stride1 , stride2 , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL ) {
    Buud_corrWW = allocate_momcorrs( stride1 , stride2 , wwNMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    error_code = FAILURE ; goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // if we are doing nonrel-chiral hadrons we switch chiral to nrel
    rotate_offdiag_2( S1 , prop1.basis , 
		      S2 , prop2.basis ) ;

    // accumulate wall sum expects both to be walls
    struct spinor SUM1 , SUM2 ;
    if( prop1.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
      sumprop( &SUM2 , S2 ) ;
    }

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;

    // strange memory access pattern threads better than what was here before
    size_t site ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      #pragma omp single nowait
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_code = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	
	size_t GSGK ;
	for( GSGK = 0 ; GSGK < stride1 ; GSGK++ ) {
	  
	  const size_t GSRC = GSGK / B_CHANNELS ;
	  const size_t GSNK = GSGK % B_CHANNELS ;

	  // zero our terms
	  size_t odc ;
	  for( odc = 0 ; odc < NSNS ; odc++ ) {
	    in[ 0 + 2 * ( odc + NSNS * GSGK ) ][ site ] = 0.0 ;
	    in[ 1 + 2 * ( odc + NSNS * GSGK ) ][ site ] = 0.0 ;
	  }
	  
	  // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
	  const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	  // precompute \gamma_t ( Cg_\mu )^{dagger} \gamma_t 
          const struct gamma Cgnu = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
          const struct gamma CgnuD = gt_Gdag_gt( Cgnu , GAMMAS[ GAMMA_T ] ) ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( in , S1[ site ] , S1[ site ] , S2[ site ] , 
				    Cgmu , CgnuD , GSGK , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop1.source == WALL ) {
	baryon_contract_walls( Buud_corrWW , SUM1 , SUM1 , SUM2 , GAMMAS , 
			       tshifted , UUD_BARYON ) ;
      }
    }

    // momentum projection 
    baryon_momentum_project( Buud_corr , in , out , forward , backward ,
			     list , NMOM , tshifted , UUD_BARYON ) ;

    // if we error we leave
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy over the propagators
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      memcpy( &S1[i] , &S1f[i] , sizeof( struct spinor ) ) ;
      memcpy( &S2[i] , &S2f[i] , sizeof( struct spinor ) ) ;
    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the baryons wall-local and maybe wall-wall
  write_baryon( Buud_corr , list , NMOM , GLU_FALSE , outfile , "uud" ) ;

  if( prop1.source == WALL ) {
    write_baryon( Buud_corrWW , wwlist , wwNMOM , GLU_TRUE , outfile , "uud" ) ;
  }

  // failure sink
 memfree :

  // free our momentum correlators
  if( NMOM != NULL ) {
    free_momcorrs( Buud_corr , stride1 , stride2 , NMOM[0] ) ;
    if( prop1.source == WALL ) {
      free_momcorrs( Buud_corrWW , stride1 , stride2 , wwNMOM[0] ) ;
    }
  }

  // free our ffts
  free_ffts( in , out , forward , backward , flat_dirac ) ;

  // free spinors
  free( S1 ) ; free( S1f ) ;
  free( S2 ) ; free( S2f ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

  // tell us how long it all took
  print_time( ) ;

  return error_code ;
}

