/**
   @file nrqcd.c
   @brief "on the fly" NRQCD code
 **/
#include "common.h"

#include "GLU_timer.h"
#include "evolve.h"

// compute our various NRQCD propagators
int
compute_nrqcd_props( struct propagator *prop ,
		     const size_t nprops )
{
  // loop propagators checking to see if any are non-rel
  // check tadpole factors are the same
  double tadref = 0.0 ;
  GLU_bool HAVE_C7C8 = GLU_FALSE ;
  GLU_bool HAVE_C11 = GLU_FALSE ;
  GLU_bool FLY_NREL = GLU_FALSE ;
  size_t i , n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].basis == NREL_CORR ) {
      FLY_NREL = GLU_TRUE ;

      // allocate the heavy propagator
      if( corr_malloc( (void**)&prop[n].H ,
		       ALIGNMENT ,
		       LVOLUME*sizeof( struct halfspinor_f ) ) != 0 ) {
	fprintf( stderr , "[NRQCD] heavy prop allocation failure\n" ) ;
	goto memfree ;
      }

      // these ones require halfspinor temporaries over LCU
      if( fabs( prop[n].NRQCD.C7 ) > NRQCD_TOL ||
	  fabs( prop[n].NRQCD.C8 ) > NRQCD_TOL ) {
	HAVE_C7C8 = GLU_TRUE ;
      }
      // this one needs two!
      if( fabs( prop[n].NRQCD.C11 ) > NRQCD_TOL ) {
	HAVE_C7C8 = GLU_TRUE ;
	HAVE_C11 = GLU_TRUE ;
      }
      
      if( fabs( tadref ) < NRQCD_TOL ) {
	tadref = prop[n].NRQCD.U0 ;
      } else {
	if( fabs( prop[n].NRQCD.U0 - tadref ) > 1E-12 ||
	    prop[n].NRQCD.U0 < 1E-12 ) {
	  fprintf( stderr , "[NRQCD] poor selection of tadpole factor "
		   "prop %zu || fac %f \n" , n , prop[n].NRQCD.U0 ) ;
	  return FAILURE ;
	}
      }
    }
  }
  // if we aren't doing any NRQCD props then we successfully do nothing
  if( FLY_NREL == GLU_FALSE ) {
    fprintf( stdout , "[NRQCD] Not computing NRQCD props on the fly\n" ) ;
    return SUCCESS ;
  }
  // check for an empty gauge field
  if( lat == NULL ) {
    fprintf( stderr , "[NRQCD] Empty gauge field, add with -c on command line\n" ) ;
    return FAILURE ;
  }  
  // otherwise we initialise all this gubbins
  struct NRQCD_fields F ;

  // initialise all temporary fields to null
  F.S = F.H = F.S1 = F.S2 = NULL ; F.Fmunu = NULL ;

  // usual allocations F.S is the prop, F.H the summed hamiltonian
  if( corr_malloc( (void**)&F.S  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ||
      corr_malloc( (void**)&F.H  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ) {
    fprintf( stderr , "[NRQCD] temporary halfspinor allocation failure\n" ) ;
    goto memfree ;
  }

  // C7 and C8 need one more temporary
  if( HAVE_C7C8 == GLU_TRUE ) {
    if( corr_malloc( (void**)&F.S1 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ) {
      fprintf( stderr , "[NRQCD] temporary halfspinor for C7 & C8 allocation failure\n" ) ;
      goto memfree ;
    }
  }
  // C11 needs another one
  if( HAVE_C11 == GLU_TRUE ) {
    if( corr_malloc( (void**)&F.S2 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ) {
      fprintf( stderr , "[NRQCD] temporary halfspinor for C11 allocation failure\n" ) ;
      goto memfree ;
    }
  }

  // allocate the gauge field
  if( corr_malloc( (void**)&F.Fmunu , ALIGNMENT , LCU*sizeof( struct field ) ) != 0 ) {
    fprintf( stderr , "[NRQCD] clover allocation failure\n" ) ;
    goto memfree ;
  }
  for( i = 0 ; i < LCU ; i++ ) {
    corr_malloc( (void**)&F.Fmunu[i].O , ALIGNMENT ,
		 ((NS-1)*(NS-2)+2*(ND-1))*sizeof( double complex* ) ) ;
    size_t j ;
    for( j = 0 ; j < (NS-1)*(NS-2)+2*(ND-1) ; j++ ) {
      corr_malloc( (void**)&F.Fmunu[i].O[j] , ALIGNMENT ,
		   NCNC*sizeof( double complex ) ) ;
    }
  }

  // initialise the timer
  start_timer() ;

#pragma omp parallel
  {
    // compute the propagators
    compute_props( prop , &F , lat , nprops , tadref ) ;
  }

  // tell us the time it took
  print_time() ;


 memfree :
  
  // memory frees of the temporary halfspinors
  if( F.S != NULL ) {
    free( F.S ) ;
  }
  if( F.S1 != NULL ) {
    free( F.S1 ) ;
  }
  if( F.S2 != NULL ) {
    free( F.S2 ) ;
  }
  if( F.H != NULL ) {
    free( F.H ) ;
  }

  // free the field strength tensor
  for( i = 0 ; i < LCU ; i++ ) {
    size_t j ;
    for( j = 0 ; j < (NS-1)*(NS-2)+2*(ND-1) ; j++ ) {
      if( F.Fmunu[i].O[j] != NULL ) {
	free( F.Fmunu[i].O[j] ) ;
      }
    }
    if( F.Fmunu[i].O != NULL ) {
      free( F.Fmunu[i].O ) ;
    }
  }
  if( F.Fmunu != NULL ) {
    free( F.Fmunu ) ;
  }
  
  return SUCCESS ;
}

int
free_nrqcd_props( struct propagator *prop ,
		  const size_t nprops )
{
  size_t n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].H != NULL ) {
      free( prop[n].H ) ;
    }
  }
  return SUCCESS ;
}
