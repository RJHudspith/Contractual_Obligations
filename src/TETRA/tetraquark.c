/**
   @file tetraquark.c
   @brief tetraquark contraction code

   This is just baryons.c for now - J ( 7 , 05 , 2015 )
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
tetraquark_diagonal( struct propagator prop ,
		     struct cut_info CUTINFO ,
		     const char *outfile )
{
  // flat dirac indices
  const int flat_dirac = B_CHANNELS * NSNS ;

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
  struct mcorr **Buud_corr = NULL , **Buuu_corr = NULL , **Buds_corr = NULL ;
  struct mcorr **Buud_corrWW = NULL , **Buuu_corrWW = NULL , **Buds_corrWW = NULL ;

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

  in = malloc( ( 2 * flat_dirac ) * sizeof( double complex ) ) ;
  int i ;
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H

  out = malloc( ( 2 * flat_dirac ) * sizeof( double complex ) ) ;
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
  Buds_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
  Buud_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
  Buuu_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop.source == WALL ) {
    Buds_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
    Buud_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
    Buuu_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop , S1 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // compute wall sum
    struct spinor SUM1 ;
    if( prop.source == WALL ) {
      sumprop( &SUM1 , S1 ) ;
    }

    // strange memory access pattern threads better than what was here before
    int site ;
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
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	
	int GSRC ;
	for( GSRC = 0 ; GSRC < ( B_CHANNELS ) ; GSRC++ ) {

	  // zero our terms
	  int odc ;
	  for( odc = 0 ; odc < NSNS ; odc++ ) {
	    in[ 0 + 2 * ( odc + NSNS * GSRC ) ][ site ] = 0.0 ;
	    in[ 1 + 2 * ( odc + NSNS * GSRC ) ][ site ] = 0.0 ;
	  }
	  
	  // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
	  const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	  // precompute \gamma_t ( Cg_\mu )^{*} \gamma_t -> \Gamma^{T} in note
	  const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( in , S1[ site ] , S1[ site ] , S1[ site ] , 
				    Cgmu , CgmuT , GSRC , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop.source == WALL ) {
	baryon_contract_walls( Buud_corrWW , Buuu_corrWW , Buds_corrWW ,
			       SUM1 , SUM1 , SUM1 , GAMMAS , t ) ;
      }
    }

    // momentum projection 
    baryon_momentum_project( Buud_corr , Buuu_corr , Buds_corr ,
			     in , out , forward , backward ,
			     list , NMOM , t ) ;

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
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the baryons wall-local and maybe wall-wall
  write_baryons( Buud_corr , Buuu_corr , Buds_corr ,
		 list , NMOM , GLU_FALSE , outfile ) ;

  if( prop.source == WALL ) {
    write_baryons( Buud_corrWW , Buuu_corrWW , Buds_corrWW ,
		   wwlist , wwNMOM , GLU_TRUE , outfile ) ;
  }

  // free our momentum correlators
  free_momcorrs( Buds_corr , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buud_corr , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buuu_corr , B_CHANNELS , NSNS , NMOM[0] ) ;

  if( prop.source == WALL ) {
    free_momcorrs( Buds_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buud_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buuu_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
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

  // free our correlators
  free_momcorrs( Buds_corr , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buud_corr , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buuu_corr , B_CHANNELS , NSNS , NMOM[0] ) ;

  if( prop.source == WALL ) {
    free_momcorrs( Buds_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buud_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buuu_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
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
