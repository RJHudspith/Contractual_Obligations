/**
   @file baryons.c
   @brief baryon contraction code
*/
#include "common.h"

#include "bar_contractions.h"  // baryon_contract_site_mom()
#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // gamma_mul_lr()
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "cut_routines.h"      // veclist
#include "gammas.h"            // make_gammas() && gamma_mmul*
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "plan_ffts.h"         // create_plans_DFT() 
#include "read_propheader.h"   // for read_propheader()
#include "spinor_ops.h"        // sumprop()

// flavour degenerate baryon contraction
int
baryons_diagonal( struct propagator prop ,
		  struct cut_info CUTINFO ,
		  const char *outfile )
{
  // flat dirac indices
  const size_t flat_dirac = ( B_CHANNELS * B_CHANNELS ) * NSNS ;

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL ;

  // momentum stuff
  int *NMOM = NULL , *wwNMOM = NULL ;
  struct veclist *list = NULL , *wwlist = NULL ;

  // fftw temporaries
  double complex **in = NULL , **out = NULL ;

  // correlators
  struct mcorr **Buuu_corr = NULL , **Buuu_corrWW = NULL ;

#ifdef HAVE_FFTW3_H
  fftw_plan *forward = NULL , *backward = NULL ;
#else
  int *forward = NULL , *backward = NULL ;
#endif

  // allocations
  if( corr_malloc( (void**)&S1 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S1f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  NMOM = malloc( sizeof( int ) ) ;
  wwNMOM = malloc( sizeof( int ) ) ;

  in = malloc( ( 2 * flat_dirac ) * sizeof( double complex* ) ) ;
  size_t i ;
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H

  out = malloc( ( 2 * flat_dirac ) * sizeof( double complex* ) ) ;
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    out[ i ] = malloc( LCU * sizeof( double complex ) ) ;
  }

  forward  = ( fftw_plan* )malloc( ( 2 * flat_dirac ) * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( ( 2 * flat_dirac ) * sizeof( fftw_plan ) ) ;

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , 2 * flat_dirac , ND-1 ) ;

  list = (struct veclist*)compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

#else

  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

#endif

  // allocate the zero veclist for our writing out purposes
  if( prop.source == WALL ) {
    wwlist = (struct veclist*)zero_veclist( wwNMOM , ND-1 , GLU_FALSE ) ;
  }

  // Define our output correlators, with B_CHANNELS channels and NSNS components
  Buuu_corr = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop.source == WALL ) {
    Buuu_corrWW = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , wwNMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  size_t t ;
  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    // shifted times
    const size_t tshifted = ( t + LT - prop.origin[ND-1] ) % LT ;

    // strange memory access pattern threads better than what was here before
    size_t site ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	
	size_t GSGK ; // combined gamma source and sink indices
	for( GSGK = 0 ; GSGK < ( B_CHANNELS * B_CHANNELS ) ; GSGK++ ) {

	  // source and sink gammas
	  const size_t GSRC = GSGK / B_CHANNELS ;
	  const size_t GSNK = GSGK % B_CHANNELS ;

	  // zero our terms
	  size_t odc ;
	  for( odc = 0 ; odc < NSNS ; odc++ ) {
	    in[ 0 + 2 * ( odc + NSNS * ( GSGK ) ) ][ site ] = 0.0 ;
	    in[ 1 + 2 * ( odc + NSNS * ( GSGK ) ) ][ site ] = 0.0 ;
	  }
	  
	  // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
	  const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	  // precompute \gamma_t ( Cg_\mu )^{dagger} \gamma_t 
	  const struct gamma Cgnu = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
	  const struct gamma CgnuD = gt_Gdag_gt( Cgnu , GAMMAS ) ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( in , S1[ site ] , S1[ site ] , S1[ site ] , 
				    Cgmu , CgnuD , GSGK , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop.source == WALL ) {
	baryon_contract_walls( Buuu_corrWW , SUM1 , SUM1 , SUM1 , GAMMAS , tshifted ,
			       UUU_BARYON ) ;
      }
    }

    // momentum projection 
    baryon_momentum_project( Buuu_corr , in , out , forward , backward ,
			     list , NMOM , tshifted , UUU_BARYON ) ;

    // if we error we leave
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
    printf("\r[BARYONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the baryons wall-local and maybe wall-wall
  write_baryon( Buuu_corr , list , NMOM , GLU_FALSE , outfile , "uuu" ) ;

  if( prop.source == WALL ) {
    write_baryon( Buuu_corrWW , wwlist , wwNMOM , GLU_TRUE , outfile , "uuu" ) ;
  }

  // free our momentum correlators
  free_momcorrs( Buuu_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;

  if( prop.source == WALL ) {
    free_momcorrs( Buuu_corrWW , B_CHANNELS * B_CHANNELS , NSNS , wwNMOM[0] ) ;
  }

  // free the "in" allocation
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    free( in[ i ] ) ;
  }
  free( in ) ;

#ifdef HAVE_FFTW3_H
  // free fftw stuff
  #pragma omp parallel for private(i)
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    fftw_destroy_plan( forward[i] ) ;
    fftw_destroy_plan( backward[i] ) ;
    fftw_free( out[i] ) ;
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; 
  fftw_cleanup( ) ; 
#endif

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free stuff
  free( S1 ) ; free( S1f ) ;

  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;

  // failure sink
 FREE_FAIL :

  // free the "in" allocation
  if( in != NULL ) {
    for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
      free( in[ i ] ) ;
    }
  }
  free( in ) ;

#ifdef HAVE_FFTW3_H
  // free fftw stuff
  if( out != NULL ) {
    for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
      free( out[ i ] ) ;
    }
  }
  if( forward != NULL ) {
    for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
      fftw_destroy_plan( forward[ i ] ) ;
    }
  }
  if( backward != NULL ) {
    for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
      fftw_destroy_plan( forward[ i ] ) ;
    }
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; 
  fftw_cleanup( ) ; 
#endif

  // free our correlators
  if( NMOM != NULL ) {
    free_momcorrs( Buuu_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    if( prop.source == WALL ) {
      free_momcorrs( Buuu_corrWW , B_CHANNELS * B_CHANNELS , NSNS , wwNMOM[0] ) ;
    }
  }

  // free spinors
  free( S1f ) ; free( S1 ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}
