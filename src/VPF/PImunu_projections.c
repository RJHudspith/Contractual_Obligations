/**
   @file PImunu_projections.c
   @brief various projections and related for the VPF
 */

#include "common.h"

#include "cut_output.h" // for the IO 

// we have the zero at index [0] as the FT is in the 0->2Pi BZ
void
subtract_zeromom( struct PIdata *__restrict data ,
		  const int NMOM[1] )
{
  printf( "[VPF] performing noisy 0 subtraction \n" ) ;
  return ;
}

// projections should go in their own file right?
void
projection( const struct PIdata *data ,
	    const double **p ,
	    const double *psq ,
	    const struct veclist *list ,
	    const int *__restrict NMOM ,
	    const char *outfile )
{
  double *trans = malloc( NMOM[0] * sizeof( double ) ) ;
  double *longitudinal = malloc( NMOM[0] * sizeof( double ) ) ;

  const double NORM = 1.0 / (double)( ND - 1 ) ;
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    const int list_idx = list[i].idx ;
    const double spsq = ( psq[i] == 0.0 ) ? 1.0 : 1.0 / psq[i] ;
    register double sumtrans = 0.0 , sumlong = 0.0 ;
    int mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	const double pmunu = p[i][mu] * p[i][nu] * spsq ;
	const double fac = ( mu != nu ) ? -pmunu : 1.0 - pmunu ;
	sumtrans += creal( data[list_idx].PI[mu][nu] ) * fac ;
	sumlong  += creal( data[list_idx].PI[mu][nu] ) * pmunu ;
      }
    }
    trans[ i ] = sumtrans * ( spsq * NORM ) ;
    longitudinal[ i ] = -sumlong * spsq ;
  }

  // write out a file
  char str[ 256 ] ;
  sprintf( str , "%s.trans.bin" , outfile ) ;
  write_momspace_data( str , NMOM , trans , list , ND ) ;

  sprintf( str , "%s.long.bin" , outfile ) ;
  write_momspace_data( str , NMOM , longitudinal , list , ND ) ;

#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    trans[ i ] += longitudinal[ i ] ; // put the 0+1 contribution in too
  }

  sprintf( str , "%s.transPlong.bin" , outfile ) ;
  write_momspace_data( str , NMOM , longitudinal , list , ND ) ;

  // free the projected data
  free( trans ) ;
  free( longitudinal ) ;

  return ;
}

// wrapper for the IO, WI checks and projections
void
momspace_data( struct PIdata *data ,
	       const double **p ,
	       const double *psq ,
	       const struct veclist *__restrict list ,
	       const int *__restrict NMOM ,
	       const char *outfile ,
	       const current_type current ,
	       const vector_axial VA ) 
{
  // perform WI correction if doing conserved currents
  switch( current ) {
  case CONSERVED_LOCAL :
    correct_WI( data , CORR_MU , list , NMOM[0] ) ;
    break ;
  case LOCAL_LOCAL :
    break ;
  }

  // tell us how much we violate the WI
  compute_WI( data , (const double**)p , list , NMOM[0] ) ;

  // and perform the projection
  char str[ 256 ] ;
  switch( VA ) {
  case VECTOR :
    switch( current ) {
    case CONSERVED_LOCAL :
      sprintf( str , "%s.CVLV" , outfile ) ;
      break ;
    case LOCAL_LOCAL :
      sprintf( str , "%s.LVLV" , outfile ) ;
      break ;
    }
    break ;
  case AXIAL :
    switch( current ) {
    case CONSERVED_LOCAL :
      sprintf( str , "%s.CALA" , outfile ) ;
      break ;
    case LOCAL_LOCAL :
      sprintf( str , "%s.LALA" , outfile ) ;
      break ;
    }
    break ;
  }

  // can we just use the other one?
  projection( data , (const double**)p , psq , list , NMOM , str ) ;

  return ;
}
