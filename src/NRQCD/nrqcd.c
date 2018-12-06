/**
   @file nrqcd.c
   @brief "on the fly" NRQCD code
 **/
#include "common.h"

#include "GLU_timer.h"
#include "evolve.h"
#include "plaqs_links.h"      // average_plaquette()

// free the allocated NRQCD fields
static int
free_NRQCD_fields( struct NRQCD_fields *F )
{
  // memory frees of the temporary halfspinors
  if( F -> S != NULL ) {
    free( F -> S ) ;
  }
  if( F -> S1 != NULL ) {
    free( F -> S1 ) ;
  }
  if( F -> S2 != NULL ) {
    free( F -> S2 ) ;
  }
  if( F -> H != NULL ) {
    free( F -> H ) ;
  }

  // free the field strength tensor
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    size_t j ;
    for( j = 0 ; j < (NS-1)*(NS-2)+2*(ND-1) ; j++ ) {
      if( F -> Fmunu[i].O[j] != NULL ) {
	free( F -> Fmunu[i].O[j] ) ;
      }
    }
    if( F -> Fmunu[i].O != NULL ) {
      free( F -> Fmunu[i].O ) ;
    }
  }
  if( F -> Fmunu != NULL ) {
    free( F -> Fmunu ) ;
  }
  return SUCCESS ;
}

static GLU_bool
is_fly_NRQCD( struct propagator *prop ,
	      double *tadref ,
	      GLU_bool *HAVE_C11 ,
	      const size_t nprops )
{
  GLU_bool FLY_NREL = GLU_FALSE ;
  size_t n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].basis == NREL_CORR ) {
      FLY_NREL = GLU_TRUE ;
      // allocate the heavy propagator
      if( prop[n].NRQCD.FWD == GLU_TRUE ) {
	if( corr_malloc( (void**)&prop[n].Hfwd ,
			 ALIGNMENT ,
			 T_NRQCD*LCU*sizeof( struct halfspinor_f ) ) != 0 ) {
	  fprintf( stderr , "[NRQCD] heavy bwd prop allocation failure\n" ) ;
	  return GLU_FALSE ;
	}
      }
      if( prop[n].NRQCD.BWD == GLU_TRUE ) {
	// allocate the heavy propagator
	if( corr_malloc( (void**)&prop[n].Hbwd ,
			 ALIGNMENT ,
			 T_NRQCD*LCU*sizeof( struct halfspinor_f ) ) != 0 ) {
	  fprintf( stderr , "[NRQCD] heavy fwd prop allocation failure\n" ) ;
	  return GLU_FALSE ;
	}
      }
      // this one needs another temporary
      if( fabs( prop[n].NRQCD.C11 ) > NRQCD_TOL ) {
	*HAVE_C11 = GLU_TRUE ;
      }
      // compare tadpole factors
      if( fabs( *tadref ) < NRQCD_TOL ) {
	*tadref = prop[n].NRQCD.U0 ;
      } else {
	if( fabs( prop[n].NRQCD.U0 - *tadref ) > 1E-12 ||
	    prop[n].NRQCD.U0 < 1E-12 ) {
	  fprintf( stderr , "[NRQCD] poor selection of tadpole factor "
		   "prop %zu || fac %f \n" , n , prop[n].NRQCD.U0 ) ;
	  return FAILURE ;
	}
      }
    }
  }
  return FLY_NREL ;
}

// compute our various NRQCD propagators
int
compute_nrqcd_props( struct propagator *prop ,
		     const size_t nprops )
{
  // check for an empty gauge field
  if( lat == NULL ) {
    fprintf( stderr , "[NRQCD] Empty gauge field, add with -c on command line\n" ) ;
    return FAILURE ;
  }
  
  // loop propagators checking to see if any are non-rel
  // check tadpole factors are the same
  double tadref = 0.0 ;
  GLU_bool HAVE_C11 = GLU_FALSE ;
  GLU_bool FLY_NREL = is_fly_NRQCD( prop , &tadref , &HAVE_C11 , nprops  ) ;
  size_t i ;

  // if we aren't doing any NRQCD props then we successfully do nothing
  if( FLY_NREL == GLU_FALSE ) {
    fprintf( stdout , "[NRQCD] Not computing NRQCD props on the fly\n" ) ;
    return SUCCESS ;
  }

  // give us a little warning about this
  if( (const size_t)T_NRQCD != LT ) {
    fprintf( stdout ,
	     "[NRQCD] evolving NRQCD to t=%zu. This differs from LT=%zu\n"
	     "[NRQCD] be careful if contracting with light props as LT gets changed to T_NRQCD\n" ,
	     (size_t)T_NRQCD , LT ) ;
  }
  
  // otherwise we initialise all this gubbins
  struct NRQCD_fields F ;

  // initialise all temporary fields to null
  F.S = F.H = F.S1 = F.S2 = NULL ; F.Fmunu = NULL ;

  // usual allocations F.S is the prop, F.H the summed hamiltonian
  if( corr_malloc( (void**)&F.S  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ||
      corr_malloc( (void**)&F.S1 , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ||
      corr_malloc( (void**)&F.H  , ALIGNMENT , LCU*sizeof( struct halfspinor ) ) != 0 ) {
    fprintf( stderr , "[NRQCD] temporary halfspinor allocation failure\n" ) ;
    goto memfree ;
  }

  // C11 needs another temporary
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

  const double plaq = av_plaquette( lat ) ;
  for( i = 0 ; i < nprops ; i++ ) {
    // set the plaquette to the computed value, disregarding the prop value
    prop[i].plaq = plaq ;
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

  // change the global temporal length to match the NRQCD one we set
  Latt.dims[ND-1] = (const size_t)T_NRQCD ;

 memfree :

  free_NRQCD_fields( &F ) ;
    
  return SUCCESS ;
}

int
free_nrqcd_props( struct propagator *prop ,
		  const size_t nprops )
{
  size_t n ;
  for( n = 0 ; n < nprops ; n++ ) {
    if( prop[n].Hfwd != NULL ) {
      free( prop[n].Hfwd ) ;
    }
    if( prop[n].Hbwd != NULL ) {
      free( prop[n].Hbwd ) ;
    }
  }
  return SUCCESS ;
}
