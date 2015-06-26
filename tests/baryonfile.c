/**
   @file baryonfile.c
   @brief little correlation file reader
   
   TODO :: momentum projections for spin 1/2, 3/2?
 */
#include "common.h"

#include "contractions.h" // dirac_trace
#include "correlators.h"  // free_momcorrs()
#include "gammas.h"       // make_gammas()
#include "reader.h"       // process file

// global lattice dimensions and stuff
struct latt_info Latt ;

// enum for the projections we allow
typedef enum { LO , L1 , L2 , L3 , L4 , L5 } bprojection ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// computes res = G * spinmatrix
static void
gamma_spinmatrix( void *res , 
		  const struct gamma G ,
		  const void *spinmatrix ) 
{
  double complex *r = (double complex*)res ;
  const double complex *d = (const double complex*)spinmatrix ;
  int i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    register double complex factor = 1.0 ;
    switch( G.g[i] ) {
    case 0 : factor = +1 ; break ;
    case 1 : factor = +I ; break ;
    case 2 : factor = -1 ; break ;
    case 3 : factor = -I ; break ;
    }
    for( j = 0 ; j < NS ; j++ ) {
      r[ j + i * NS ] = factor * d[ j + G.ig[i] * NS ] ;
    }
  }
}

// atomically add two spinmatrices
static void
atomic_add_spinors( void *res ,
		    const void *D )
{
  double complex *r = (double complex*)res ;
  const double complex *d = (const double complex*)D ;
  int i ;
  for( i = 0 ; i < NSNS ; i++ ) {
    r[ i ] += d[ i ] ;
  }
}

// spin 1/2 projection? summed into D
void
spinhalf_project( double complex *D , 
		  const struct mcorr **corr ,
		  const struct gamma *GAMMA ,
		  const size_t i ,
		  const size_t p ,
		  const size_t t )
{
  double complex *tmp1 = NULL , *tmp2 = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&tmp2 , 16 , NSNS * sizeof( double complex ) ) ;

  int j , d ;
  for( d = 0 ; d < NSNS ; d++ ) {
    D[ d ] = 0.0 ;
  }

  // perform 1/3 \gamma_i \gamma_j O_j
  for( j = 0 ; j < ND-1 ; j++ ) {
    for( d = 0 ; d < NSNS ; d++ ) {
      tmp1[ d ] = corr[ j+i*B_CHANNELS ][ d ].mom[ p ].C[ t ] ;
    }
    struct gamma gij ; gamma_mmul( &gij , GAMMA[i] , GAMMA[j] ) ;
    gamma_spinmatrix( tmp2 , gij , tmp1 ) ;
    atomic_add_spinors( D , tmp2 ) ;
  }

  // normalise
  for( d = 0 ; d < NSNS ; d++ ) {
    D[ d ] /= 3.0 ;
  }
  free( tmp1 ) ;
  free( tmp2 ) ;
  return ;
}

// spin 1/2 projection? summed into D
void
spinthreehalf_project( double complex *D , 
		       const struct mcorr **corr ,
		       const struct gamma *GAMMA ,
		       const size_t i ,
		       const size_t p ,
		       const size_t t )
{
  double complex *tmp1 = NULL , *tmp2 = NULL ;
  corr_malloc( (void**)&tmp1 , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&tmp2 , 16 , NSNS * sizeof( double complex ) ) ;

  int j , d ;
  for( d = 0 ; d < NSNS ; d++ ) {
    D[ d ] = 0.0 ;
  }

  // perform delta_{ij} - 1/3 \gamma_i \gamma_j O_j
  for( j = 0 ; j < ND-1 ; j++ ) {
    for( d = 0 ; d < NSNS ; d++ ) {
      tmp1[ d ] = corr[ j+i*B_CHANNELS ][ d ].mom[ p ].C[ t ] ;
    }
    if( j == i ) {
      atomic_add_spinors( D , tmp1 ) ;
    }
    struct gamma gij ; gamma_mmul( &gij , GAMMA[i] , GAMMA[j] ) ;
    gamma_spinmatrix( tmp2 , gij , tmp1 ) ;
    // normalise
    for( d = 0 ; d < NSNS ; d++ ) {
      tmp2[ d ] /= -3.0 ;
    }
    atomic_add_spinors( D , tmp2 ) ;
  }

  free( tmp1 ) ;
  free( tmp2 ) ;
  return ;
}

// returns the correlator at some momentum
double complex*
baryon_project( const struct mcorr **corr ,
		const struct gamma *GAMMA ,
		const size_t GSRC ,
		const size_t GSNK ,
		const size_t p ,
		const bprojection projection ) 
{
  // D is a temporary spinmatrix
  double complex *D = NULL , *result = NULL ;

  corr_malloc( (void**)&D , 16 , NSNS * sizeof( double complex ) ) ;
  corr_malloc( (void**)&result , 16 , L0 * sizeof( double complex ) ) ;

  // build these just in case
  struct gamma g0g5 ; gamma_mmul( &g0g5 , GAMMA[ 0 ] , GAMMA[ 5 ] ) ;
  struct gamma g3g0g5 ; gamma_mmul( &g3g0g5 , GAMMA[ 3 ] , g0g5 ) ;

  // loop times
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {
    
    // set D
    int i ;
    for( i = 0 ; i < NSNS ; i++ ) {
      D[ i ] = corr[ GSNK + B_CHANNELS * GSRC ][ i ].mom[ p ].C[ t ] ;
    }

    switch( projection ) {
    case LO : // does 1/4( g_I + g_3 - i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) + 
			     dirac_trace( GAMMA[ GAMMA_3 ] , D ) - 
			     I * dirac_trace( g0g5 , D ) - 
			     I * dirac_trace( g3g0g5 , D ) ) ;
      break ;
    case L1 : // does 1/4( g_I - g_3 + i*g_0g_5 - i*g_3g_0g_5
      result[ t ] = 0.25 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) - 
			     dirac_trace( GAMMA[ GAMMA_3 ] , D ) +
			     I * dirac_trace( g0g5 , D ) - 
			     I * dirac_trace( g3g0g5 , D ) ) ;
      break ;
    case L2 : // does 1/4( g_I - g_3 - i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) - 
			     dirac_trace( GAMMA[ GAMMA_3 ] , D ) -
			     I * dirac_trace( g0g5 , D ) + 
			     I * dirac_trace( g3g0g5 , D ) ) ;
      break ;
    case L3 : // does 1/4( g_I + g_3 + i*g_0g_5 + i*g_3g_0g_5
      result[ t ] = 0.25 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) + 
			     dirac_trace( GAMMA[ GAMMA_3 ] , D ) +
			     I * dirac_trace( g0g5 , D ) + 
			     I * dirac_trace( g3g0g5 , D ) ) ;
      break ;
    case L4 :
      result[ t ] = 0.5 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) +
			    dirac_trace( GAMMA[ GAMMA_3 ] , D ) ) ;
      break ;
    case L5 :
      result[ t ] = 0.5 * ( dirac_trace( GAMMA[ IDENTITY ] , D ) - 
			    dirac_trace( GAMMA[ GAMMA_3 ] , D ) ) ;
      break ;
    }
  }

  free( D ) ;
  return result ;
}

// little code for accessing elements of our baryon correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 4 ) {
    return printf( "usage ./BARYONS {correlator file} BASIS PROJECTION GSRC,GSNK,px,py,pz ... \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return FAILURE ;
  }
  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // structs
  struct gamma *GAMMAS = NULL ;
  struct mcorr **corr = NULL ;
  int **momentum = NULL ;

  // set the basis
  proptype basis = CHIRAL ;
  if( are_equal( argv[2] , "NREL" ) ) {
    basis = NREL ;
  } else if( are_equal( argv[2] , "STATIC" ) ) {
    basis = STATIC ;
  } else if( are_equal( argv[2] , "CHIRAL" ) ) {
    basis = CHIRAL ;
  } else {
    printf( "[INPUTS] I don't understand your basis %s \n" , argv[2] ) ;
    goto memfree ;
  }

  // set the projection
  bprojection projection = LO ;
  if( are_equal( argv[3] , "L0" ) ) {
    projection = LO ;
  } else if( are_equal( argv[3] , "L1" ) ) {
    projection = L1 ;
  } else if( are_equal( argv[3] , "L2" ) ) {
    projection = L2 ;
  } else if( are_equal( argv[3] , "L3" ) ) {
    projection = L3 ;
  } else if( are_equal( argv[3] , "L4" ) ) {
    projection = L4 ;
  } else if( are_equal( argv[3] , "L5" ) ) {
    projection = L5 ;
  } else {
    printf( "[INPUTS] I don't understand your projection %s \n" , argv[3] ) ;
    goto memfree ;
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , basis ) == FAILURE ) {
    goto memfree ;
  }

  // number of correlators printed to stdout
  int corrs_written = 0 ;

  // is defined in reader.c
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;

  if( corr == NULL ) goto memfree ;

  // loop the ones we want
  int i ;
  for( i = 4 ; i < ( argc ) ; i++ ) {
    // tokenize argv into the correlators people want
    char *tok1 = strtok( (char*)argv[i] , "," ) ;
    if( tok1 == NULL ) break ;
    const int idx1 = atoi( tok1 ) ;
    if( idx1 >= B_CHANNELS || idx1 < 0 ) { 
      printf( "[Momcorr] Non-sensical source index %d \n" , idx1 ) ;
      break ;
    } 
    char *tok2 = strtok( NULL , "," ) ;
    if( tok2 == NULL ) break ;
    const int idx2 = atoi( tok2 ) ;
    if( idx2 >= B_CHANNELS || idx2 < 0 ) { 
      printf( "[Momcorr] Non-sensical sink index %d \n" , idx2 ) ;
      break ;
    } 

    // initialise to 0
    int moms[ ND - 1 ] , mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      moms[ mu ] = 0 ;
    }

    printf( "[Momcorr] searching for momenta (" ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      char *ptok = strtok( NULL , "," ) ;
      if( ptok == NULL ) break ;
      moms[ mu ] = (int)atoi( ptok ) ;
      printf( " %d " , moms[ mu ] ) ;
    }
    printf( ") \n" ) ;

    // find the correlator in the list
    const int matchmom = find_desired_mom( (const int**)momentum , moms , 
					   (int)NMOM[0] ) ;
    if( matchmom == FAILURE ) {
      printf( "[Momcorr] Unable to find desired momentum ... Leaving \n" ) ;
      break ;
    }

    printf( "[Momcorr] match ( %d %d %d ) \n" , momentum[ matchmom ][ 0 ] ,
	    momentum[ matchmom ][ 1 ] ,  momentum[ matchmom ][ 2 ] ) ;

    printf( "[Momcorr] Correlator [ Source :: %d | Sink :: %d ] \n\n" , 
	    idx1 , idx2 ) ;

    // do the projection
    const double complex *C = baryon_project( (const struct mcorr**)corr , 
					      GAMMAS ,
					      idx1 , idx2 , matchmom ,
					      projection ) ;
    
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      printf( "CORR %d %1.12e %1.12e\n" , t ,
	      creal( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ,
	      cimag( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ) ;
    }

    free( (void*)C ) ;
    corrs_written++ ;
    //
    printf( "\n" ) ;
  }

  // if we don't have a match or didn't specify gammas give the momentum
  // list as an option
  if( corrs_written == 0 ) {
    write_momlist( (const int **)momentum , NMOM[ 0 ] ) ;
  }

 memfree :

  // free the gamma matrices
  free( GAMMAS ) ;

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  size_t p ;
  for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
    free( momentum[ p ] ) ;
  }
  free( momentum ) ;

  fclose( infile ) ;

  return 0 ;
}
