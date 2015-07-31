/**
   @file baryons2.c
   @brief two flavour baryon contraction code
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
baryons_2fdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    struct cut_info CUTINFO ,
		    const char *outfile )
{
  // flat dirac indices
  const int flat_dirac = ( B_CHANNELS * B_CHANNELS ) * NSNS ;

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
  if( corr_malloc( (void**)&S2 , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
    goto FREE_FAIL ;
  }
  if( corr_malloc( (void**)&S2f , 16 , VOL3 * sizeof( struct spinor ) ) != 0 ) {
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

  NMOM = malloc( sizeof( int ) ) ;
  wwNMOM = malloc( sizeof( int ) ) ;

  in = malloc( ( 2 * flat_dirac ) * sizeof( double complex* ) ) ;
  int i ;
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    in[ i ] = calloc( LCU , sizeof( double complex ) ) ;
  }

#ifdef HAVE_FFTW3_H

  out = malloc( ( 2 * flat_dirac ) * sizeof( double complex* ) ) ;
  for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
    out[ i ] = fftw_malloc( LCU * sizeof( double complex ) ) ;
  }

  forward = ( fftw_plan* )malloc( ( 2 * flat_dirac ) * sizeof( fftw_plan ) ) ; 
  backward = ( fftw_plan* )malloc( ( 2 * flat_dirac ) * sizeof( fftw_plan ) ) ;

  // create spatial volume fftw plans
  create_plans_DFT( forward , backward , in , out , 2 * flat_dirac , ND-1 ) ;

  list = (struct veclist*)compute_veclist( NMOM , CUTINFO , ND-1 , GLU_FALSE ) ;

#else

  list = (struct veclist*)zero_veclist( NMOM , ND-1 , GLU_FALSE ) ;

#endif

  // allocate the zero veclist for our writing out purposes
  if( prop1.source == WALL ) {
    wwlist = (struct veclist*)zero_veclist( wwNMOM , ND-1 , GLU_FALSE ) ;
  }

  // Define our output correlators, with B_CHANNELS channels and NSNS components
  Buds_corr = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
  Buud_corr = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
  Buuu_corr = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL ) {
    Buds_corrWW = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    Buud_corrWW = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    Buuu_corrWW = allocate_momcorrs( B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
  }

  // read in the first timeslice
  if( read_prop( prop1 , S1 ) == FAILURE ||
      read_prop( prop2 , S2 ) == FAILURE ) {
    goto FREE_FAIL ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

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

    // multiple time source support
    const size_t tshifted = ( t + LT - prop1.origin[ND-1] ) % LT ;

    // strange memory access pattern threads better than what was here before
    int site ;
    int error_flag = SUCCESS ;
    #pragma omp parallel
    {
      #pragma omp master
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop1 , S1f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      #pragma omp single nowait
      {
	if( t < ( LT - 1 ) ) {
	  if( read_prop( prop2 , S2f ) == FAILURE ) {
	    error_flag = FAILURE ;
	  }
	}
      }
      // Loop over spatial volume threads better
      #pragma omp for private(site) schedule(dynamic)
      for( site = 0 ; site < LCU ; site++ ) {
	
	int GSGK ;
	for( GSGK = 0 ; GSGK < ( B_CHANNELS * B_CHANNELS ) ; GSGK++ ) {
	  
	  const int GSRC = GSGK / B_CHANNELS ;
	  const int GSNK = GSGK % B_CHANNELS ;

	  // zero our terms
	  int odc ;
	  for( odc = 0 ; odc < NSNS ; odc++ ) {
	    in[ 0 + 2 * ( odc + NSNS * GSGK ) ][ site ] = 0.0 ;
	    in[ 1 + 2 * ( odc + NSNS * GSGK ) ][ site ] = 0.0 ;
	  }
	  
	  // precompute Cg_\mu is the product, gamma_t gamma_y gamma_[GSRC]
	  const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	  // precompute \gamma_t ( Cg_\nu )^{*} \gamma_t -> \Gamma^{T} in note
	  const struct gamma Cgnu = CGmu( GAMMAS[ GSNK ] , GAMMAS ) ;
	  const struct gamma CgnuT = CGmuT( Cgnu , GAMMAS ) ;
	  
	  // Wall-Local
	  baryon_contract_site_mom( in , S1[ site ] , S1[ site ] , S2[ site ] , 
				    Cgmu , CgnuT , GSGK , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop1.source == WALL ) {
	baryon_contract_walls( Buud_corrWW , Buuu_corrWW , Buds_corrWW ,
			       SUM1 , SUM1 , SUM2 , GAMMAS , tshifted ) ;
      }
    }

    // momentum projection 
    baryon_momentum_project( Buud_corr , Buuu_corr , Buds_corr ,
			     in , out , forward , backward ,
			     list , NMOM , tshifted ) ;

    // if we error we leave
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
    printf("\r[BARYONS] done %.f %%", (t+1)/((LT)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the baryons wall-local and maybe wall-wall
  write_baryons( Buud_corr , Buuu_corr , Buds_corr ,
		 list , NMOM , GLU_FALSE , outfile ) ;

  if( prop1.source == WALL ) {
    write_baryons( Buud_corrWW , Buuu_corrWW , Buds_corrWW ,
		   wwlist , wwNMOM , GLU_TRUE , outfile ) ;
  }

  // free our momentum correlators
  free_momcorrs( Buds_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buud_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buuu_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;

  if( prop1.source == WALL ) {
    free_momcorrs( Buds_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buud_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buuu_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
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
  free( S2 ) ; free( S2f ) ;

  // free the gammas
  free( GAMMAS ) ;

  // rewind file and read header again
  rewind( prop1.file ) ; read_propheader( &prop1 ) ;
  rewind( prop2.file ) ; read_propheader( &prop2 ) ;

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
      fftw_destroy_plan( forward[i] ) ;
    }
  }
  if( backward != NULL ) {
    for( i = 0 ; i < ( 2 * flat_dirac ) ; i++ ) {
      fftw_destroy_plan( backward[i] ) ;
    }
  }
  free( forward )  ; free( backward ) ; 
  fftw_free( out ) ; 
  fftw_cleanup( ) ; 
#endif

  // free our momentum correlators
  if( NMOM != NULL ) {
    free_momcorrs( Buds_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buud_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    free_momcorrs( Buuu_corr , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    if( prop1.source == WALL ) {
      free_momcorrs( Buds_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
      free_momcorrs( Buud_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
      free_momcorrs( Buuu_corrWW , B_CHANNELS * B_CHANNELS , NSNS , NMOM[0] ) ;
    }
  }

  // free spinors
  free( S1 ) ; free( S1f ) ;
  free( S2 ) ; free( S2f ) ;

  // free momentum stuff
  free( NMOM ) ; free( (void*)list ) ;
  free( wwNMOM ) ; free( (void*)wwlist ) ;

  // free the gammas
  free( GAMMAS ) ;

  return FAILURE ;
}

