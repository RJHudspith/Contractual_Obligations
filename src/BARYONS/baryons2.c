/**
   @file baryons.c
   @brief baryon contraction code

   TODO :: flatten some part of the src/snk indices and master/slave the IO
   seems tough will need to think about it ...
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
  const int flat_dirac = B_CHANNELS * NSNS ;

  // gamma matrices
  struct gamma *GAMMAS = NULL ;

  // and our spinors
  struct spinor *S1 = NULL , *S1f = NULL , *S2 = NULL , *S2f = NULL ;

  // momentum stuff
  int *NMOM = NULL , *wwNMOM = NULL ;
  struct veclist *list = NULL , *wwlist = NULL ;

  double complex **in = NULL ;

  // correlators
  struct mcorr **Buud_corr = NULL , **Buuu_corr = NULL , **Buds_corr = NULL ;
  struct mcorr **Buud_corrWW = NULL , **Buuu_corrWW = NULL , **Buds_corrWW = NULL ;

#ifdef HAVE_FFTW3_H
  // allocate these, they are getting FFTD, or summed
  double complex **out = NULL ;
  // fftw plans
  fftw_plan *forward = NULL , *backward = NULL ; 
#endif

  // allocations
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
  Buds_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
  Buud_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
  Buuu_corr = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;

  // allocate the walls if we are using wall source propagators
  if( prop1.source == WALL ) {
    Buds_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
    Buud_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
    Buuu_corrWW = allocate_momcorrs( B_CHANNELS , NSNS , NMOM[0] ) ;
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

    // strange memory access pattern threads better than what was here before
    int site ;
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
	  baryon_contract_site_mom( in , S1[ site ] , S1[ site ] , S2[ site ] , 
				    Cgmu , CgmuT , GSRC , site ) ;
	}
      }
      // loop over open indices performing wall contraction
      if( prop1.source == WALL ) {
	// accumulate the sums with open dirac indices
	int GSRC ;
        #pragma omp parallel for private(GSRC) schedule(dynamic)
	for( GSRC = 0 ; GSRC < ( B_CHANNELS ) ; GSRC++ ) {
	  // allocate these locally
	  double complex **term = malloc( 2 * sizeof( double complex* ) ) ;
	  term[0] = calloc( NSNS , sizeof( double complex ) ) ;
	  term[1] = calloc( NSNS , sizeof( double complex ) ) ;
	  // recompute teh Cgmus, these are basically free to calculate
	  const struct gamma Cgmu = CGmu( GAMMAS[ GSRC ] , GAMMAS ) ;
	  const struct gamma CgmuT = CGmuT( Cgmu , GAMMAS ) ;
	  baryon_contract_site( term , SUM1 , SUM1 , SUM2 , Cgmu , CgmuT ) ;
	  // wall contractions project to zero spatial momentum explicitly
	  int odc ;
	  for( odc = 0 ; odc < NSNS ; odc++ ) {
	    Buds_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = term[0][odc] ;
	    Buud_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = term[0][odc] + term[1][odc] ;
	    Buuu_corrWW[ GSRC ][ odc ].mom[ 0 ].C[ t ] = 2 * term[0][odc] + 4 * term[1][odc] ;
	  }
	  free( term[0] ) ; free( term[1] ) ; free( term ) ;
	}
      }
    }

    // perform the FFTS separately here
    int GSodc ; // flatteded open dirac indices
    #pragma omp parallel for private(GSodc)
    for( GSodc = 0 ; GSodc < ( flat_dirac ) ; GSodc++ ) {
      const int GSRC = GSodc / NSNS ;
      const int odc = GSodc % NSNS ;
      const int idx = 2 * GSodc ;
      #ifdef HAVE_FFTW3_H
      fftw_execute( forward[ 0 + idx ] ) ;
      fftw_execute( forward[ 1 + idx ] ) ;
      const double complex *sum1 = out[ 0 + idx ] ;
      const double complex *sum2 = out[ 1 + idx ] ;
      int p ;
      for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
	const int lid = list[ p ].idx ;
	Buds_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = sum1[ lid ] ;
	Buud_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = sum1[ lid ] + sum2[ lid ] ;
	Buuu_corr[ GSRC ][ odc ].mom[ p ].C[ t ] = 2 * sum1[ lid ] + 4 * sum2[ lid ] ;
      }
      #else
      register double complex sum1 = 0.0 , sum2 = 0.0 ;
      int site ;
      for( site = 0 ; site < LCU ; site++ ) {
	sum1 += in[ 0 + idx ][ site ] ;
	sum2 += in[ 1 + idx ][ site ] ;
      }
      Buds_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = sum1 ;
      Buud_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = sum1 + sum2 ;
      Buuu_corr[ GSRC ][ odc ].mom[ 0 ].C[ t ] = 2 * sum1 + 4 * sum2 ;
      #endif
    }

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
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s.uds" , outfile ) ;
  write_momcorr( outstr , (const struct mcorr**)Buds_corr , 
		 list , B_CHANNELS , NSNS , NMOM ) ;
  free_momcorrs( Buds_corr , B_CHANNELS , NSNS , NMOM[ 0 ] ) ;

  sprintf( outstr , "%s.uud" , outfile ) ;
  write_momcorr( outstr , (const struct mcorr**)Buud_corr , 
		 list , B_CHANNELS , NSNS , NMOM ) ;
  free_momcorrs( Buud_corr , B_CHANNELS , NSNS , NMOM[ 0 ] ) ;

  sprintf( outstr , "%s.uuu" , outfile ) ;
  write_momcorr( outstr , (const struct mcorr**)Buuu_corr , 
		 list , B_CHANNELS , NSNS , NMOM ) ;
  free_momcorrs( Buuu_corr , B_CHANNELS , NSNS , NMOM[ 0 ] ) ;

  // and the walls
  if( prop1.source == WALL ) {
    sprintf( outstr , "%s.uds.WW" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)Buds_corrWW , 
		   wwlist , B_CHANNELS , NSNS , wwNMOM ) ;
    free_momcorrs( Buds_corrWW , B_CHANNELS , NSNS , wwNMOM[ 0 ] ) ;
    sprintf( outstr , "%s.uud.WW" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)Buud_corrWW , 
		   wwlist , B_CHANNELS , NSNS , wwNMOM ) ;
    free_momcorrs( Buud_corrWW , B_CHANNELS , NSNS , wwNMOM[ 0 ] ) ;
    sprintf( outstr , "%s.uuu.WW" , outfile ) ;
    write_momcorr( outstr , (const struct mcorr**)Buuu_corrWW , 
		   wwlist , B_CHANNELS , NSNS , wwNMOM ) ;
    free_momcorrs( Buuu_corrWW , B_CHANNELS , NSNS , wwNMOM[ 0 ] ) ;
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
  free_momcorrs( Buds_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buud_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;
  free_momcorrs( Buuu_corrWW , B_CHANNELS , NSNS , NMOM[0] ) ;

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

