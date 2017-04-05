/**
   @file mesons.c
   @brief dispersion relation computation for mesons
 */
#include "common.h"

#include "basis_conversions.h" // chiral->nrel
#include "contractions.h"      // meson contract
#include "correlators.h"       // for allocate_corrs and free_corrs
#include "gammas.h"            // gt_Gdag_gt()
#include "io.h"                // read_prop
#include "progress_bar.h"      // progress_bar()
#include "setup.h"             // free_ffts() ..
#include "spinor_ops.h"        // sumprop()
#include "spinmatrix_ops.h"    // gammaspinmatrix_trace()

// number of propagators
#define Nprops (1)

// computes flavour-diagonal correlators
int
mesons_diagonal( struct propagator prop1 ,
		 const struct cut_info CUTINFO ,
		 const char *outfile )
{
  // counters
  const size_t stride1 = NSNS ;
  const size_t stride2 = NSNS ;

  // flat dirac indices are all colors and all single gamma combinations
  const size_t flat_dirac = stride1 * stride2 ;

  // error flag for if the code messes up
  int error_code = SUCCESS ;

  // loop counters
  size_t t , GSGK ;

  // initialise our measurement struct
  struct propagator prop[ Nprops ] = { prop1 } ;
  struct measurements M ;
  if( init_measurements( &M , prop , Nprops , CUTINFO ,
			 stride1 , stride2 , flat_dirac ) == FAILURE ) {
    fprintf( stderr , "[MESONS] failure to initialise measurements\n" ) ;
    error_code = FAILURE ; goto memfree ;
  }

  // initially read in a timeslice
  #pragma omp parallel
  {
    read_ahead( prop , M.S , &error_code , Nprops ) ;
  }
  if( error_code == FAILURE ) {
    goto memfree ;
  }

  // Time slice loop 
  for( t = 0 ; t < LT ; t++ ) {

    // compute wall-wall sum
    if( M.is_wall == GLU_TRUE ) {
      sumwalls( M.SUM , (const struct spinor**)M.S , Nprops ) ;
    }
    
    // for multiple time sources
    const size_t tshifted = ( t - prop1.origin[ ND-1 ] + LT ) % LT ;

    // master-slave the IO and perform each FFT (if available) in parallel
    #pragma omp parallel
    {
      if( t < ( LT - 1 ) ) {
	error_code = read_ahead( prop , M.Sf , &error_code , Nprops ) ;
      }
      /*
      // parallelise the furthest out loop :: flatten the gammas
      #pragma omp for private(GSGK) schedule(dynamic)
      for( GSGK = 0 ; GSGK < flat_dirac ; GSGK++ ) {
	const size_t GSRC = GSGK / stride1 ;
	const size_t GSNK = GSGK % stride2 ;
	const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSNK ] , 
						       M.GAMMAS[ GAMMA_T ] ) ;
	// loop spatial hypercube
	size_t site ;
	for( site = 0 ; site < LCU ; site++ ) {
	  M.in[ GSGK ][ site ] = 
	    meson_contract( gt_GSNKdag_gt  , M.S[0][ site ] , 
			    M.GAMMAS[ GSRC ] , M.S[0][ site ] ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
	// correlator computed just out of the summed walls
	if( M.is_wall == GLU_TRUE ) {
	  M.wwcorr[ GSRC ][ GSNK ].mom[0].C[ tshifted ] =	\
	    meson_contract( gt_GSNKdag_gt  , M.SUM[0] , 
			    M.GAMMAS[ GSRC ] , M.SUM[0] ,
			    M.GAMMAS[ GAMMA_5 ] ) ;
	}
      }
      */
      // if t == source then we compute disconnected loops
      if( ( t == prop1.origin[ ND-1 ] ) && ( M.is_wall == GLU_TRUE ) ) {

	size_t i , d , c ;
	register double complex C = 0.0 ;
	for( i = 0 ; i < 1 ; i++ ) {
	  struct spinor Adj ;
	  full_adj( &Adj , M.S[0][ i ] , M.GAMMAS[ GAMMA_5 ] ) ; // adjoint
	  gamma_mul_r( &Adj , M.GAMMAS[ GAMMA_5 ] ) ;
	  for( d = 0 ; d < NS ; d++ ) {
	    for( c = 0 ; c < NC ; c++ ) {
	      C += Adj.D[d][d].C[c][c] ;
	    }
	  }
	}
	printf( "%e %e \n" , creal( C ) , cimag( C ) ) ;


	struct spinor AdjSum ;
	double complex *S ;
	corr_malloc( (void**)&S , 16 , NSNS * sizeof( double complex ) ) ;
	full_adj( &AdjSum , M.SUM[0] , M.GAMMAS[ GAMMA_5 ] ) ; // adjoint
	colortrace_spinor( S , &AdjSum ) ; // trace out color indices
        //#pragma omp for private(GSGK) schedule(dynamic)
	//for( GSGK = 0 ; GSGK < stride1 ; GSGK++ ) {
	GSGK = 5 ;
	  const struct gamma gt_GSNKdag_gt = gt_Gdag_gt( M.GAMMAS[ GSGK ] , 
							 M.GAMMAS[ GAMMA_T ] ) ;
	  double complex C1 = gammaspinmatrix_trace( M.GAMMAS[ GSGK ] , S ) ;

	  printf( "%zu %1.12e %1.12e \n" , GSGK , creal( C1 ) , cimag( C1 ) ) ; 
	  
	  double complex C2 = gammaspinmatrix_trace( gt_GSNKdag_gt , S ) ;

	  printf( "%zu %1.12e %1.12e \n" , GSGK , creal( C2 ) , cimag( C2 ) ) ; 

	  C1*=C2 ;
	  printf( "%f %f \n\n" , creal( C1 ) , cimag( C1 ) ) ;
	  //}
	printf( "%zu \n" , LVOLUME ) ;
	free( S ) ;

	exit(1) ;
      }
      // end of parallel loop
    }

    // compute the contracted correlator
    compute_correlator( &M , stride1 , stride2 , tshifted ,
			CUTINFO.configspace ) ;

    // to err is human
    if( error_code == FAILURE ) {
      goto memfree ;
    }

    // copy Sf into S
    copy_props( &M , Nprops ) ;

    // status of the computation
    progress_bar( t , LT ) ;
  }

  // write out the ND-1 momentum-injected correlator and maybe the wall
  write_momcorr( outfile , (const struct mcorr**)M.corr , 
		 M.list , stride1 , stride2 , M.nmom , "" ) ;
  if( M.is_wall == GLU_TRUE ) {
    write_momcorr( outfile , (const struct mcorr**)M.wwcorr ,
		   M.wwlist , stride1 , stride2 , M.wwnmom , "ww" ) ;
  }

 memfree :

  // free our measurement struct
  free_measurements( &M , Nprops , stride1 , stride2 , flat_dirac ) ;

  return error_code ;
}

// clean up number of props
#undef Nprops
