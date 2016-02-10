/**
   @file baryonfile.c
   @brief little correlation file reader
   
   TODO :: momentum projections for spin 1/2, 3/2?
 */
#include "common.h"

#include "bar_projections.h"  // enums and projections
#include "correlators.h"      // free_momcorrs()
#include "gammas.h"           // make_gammas()
#include "reader.h"           // process file

// global lattice dimensions and stuff
struct latt_info Latt ;

// usage the code expects
static int
usage( void )
{
  return printf( "[BARYON] usage :: ./BARYONS {correlator file} "
		 " BASIS SPIN PARITY TFLIP Lx,Ly,Lz,Lt.. {outfile} \n" ) ;
}

// provide a help function
static int
help( void ) 
{
  printf( "\n\n" ) ;
  usage( ) ;
  printf( "\nI list the various options for this code, I give the \n" 
	  "OPTION :: {then a brief description} 'AND'|'THEN'|'THE'|'POSSIBLE'|'OPTIONS' \n\n" ) ;
  printf( "BASIS :: {gamma basis} 'CHIRAL'|'STATIC'|'NONREL' \n" ) ;
  printf( "SPIN :: {spin projection} 'NONE'|'1/2_11'|'1/2_12'|'1/2_21'|'1/2_22'|'3/2' \n" ) ;
  printf( "PARITY :: {parity projection} 'L0'|'L1'|'L2'|'L3'|'L4'|'L5' \n" ) ;
  printf( "TFLIP :: {time axis flip} true|false" ) ;
  return printf( "\n" ) ;
}

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// enum for the options
typedef enum {
  INFILE = 1 , GAMMA_BASIS = 2 , SPIN_PROJ = 3 , PARITY_PROJ = 4 , TIME_FLIP = 5 , DIMENSIONS = 6 , OUTFILE = 7 
} command_line_options ;

// little code for accessing elements of our baryon correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  // have a help function now to avoid confusion
  if( argc != 8 ) {
    if( argc == 1 ) {
      return usage( ) ;
    } else if( are_equal( argv[ INFILE ] , "--help" ) ) {
      return help( ) ;
    } else {
      return usage( ) ;
    }
  }

  // read the correlation file
  FILE *infile = fopen( argv[ INFILE ] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[ INFILE ] ) ;
    return FAILURE ;
  }
  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  struct gamma *GAMMAS = NULL ;     // Gamma matrices
  struct mcorr **corr = NULL ;      // read-in correlator
  struct veclist *momentum = NULL ; // momentum list
  struct mcorr **proj_corr = NULL ; // projected correlator

  // set the basis
  proptype basis = CHIRAL ;
  if( are_equal( argv[ GAMMA_BASIS ] , "NREL" ) ) {
    basis = NREL ;
  } else if( are_equal( argv[ GAMMA_BASIS ] , "STATIC" ) ) {
    basis = STATIC ;
  } else if( are_equal( argv[ GAMMA_BASIS ] , "CHIRAL" ) ) {
    basis = CHIRAL ;
  } else {
    printf( "[INPUTS] I don't understand your basis %s \n" , 
	    argv[ GAMMA_BASIS ] ) ;
    goto memfree ;
  }

  // set the spin projection
  spinhalf spin_proj = NONE ;
  if( are_equal( argv[ SPIN_PROJ ] , "1/2_11" ) ) {
    printf( "[SPIN] Performing the 11 spin-1/2 projection\n" ) ;
    spin_proj = OneHalf_11 ;
  } else if( are_equal( argv[ SPIN_PROJ ] , "1/2_12" ) ) {
    printf( "[SPIN] Performing the 12 spin-1/2 projection\n" ) ;
    spin_proj = OneHalf_12 ;
  } else if( are_equal( argv[ SPIN_PROJ ] , "1/2_21" ) ) {
    printf( "[SPIN] Performing the 21 spin-1/2 projection\n" ) ;
    spin_proj = OneHalf_21 ;
  } else if( are_equal( argv[ SPIN_PROJ ] , "1/2_22" ) ) {
    printf( "[SPIN] Performing the 22 spin-1/2 projection\n" ) ;
    spin_proj = OneHalf_22 ;
  } else if( are_equal( argv[ SPIN_PROJ ] , "3/2" ) ) {
    printf( "[SPIN] Performing the spin-3/2 projection\n" ) ;
    spin_proj = ThreeHalf ;
  } else {
    printf( "[SPIN] NOT performing a spin-projection\n" ) ;
  }

  // set the parity projection
  bprojection parity_proj = L0 ;
  if( are_equal( argv[ PARITY_PROJ ] , "L0" ) ) {
    printf( "[PARITY] Performing an L0 projection\n" ) ;
    parity_proj = L0 ;
  } else if( are_equal( argv[ PARITY_PROJ ] , "L1" ) ) {
    printf( "[PARITY] Performing an L1 projection\n" ) ;
    parity_proj = L1 ;
  } else if( are_equal( argv[ PARITY_PROJ ] , "L2" ) ) {
    printf( "[PARITY] Performing an L2 projection\n" ) ;
    parity_proj = L2 ;
  } else if( are_equal( argv[ PARITY_PROJ ] , "L3" ) ) {
    printf( "[PARITY] Performing an L3 projection\n" ) ;
    parity_proj = L3 ;
  } else if( are_equal( argv[ PARITY_PROJ ] , "L4" ) ) {
    printf( "[PARITY] Performing an L4 projection\n" ) ;
    parity_proj = L4 ;
  } else if( are_equal( argv[ PARITY_PROJ ] , "L5" ) ) {
    printf( "[PARITY] Performing an L5 projection\n" ) ;
    parity_proj = L5 ;
  } else {
    printf( "[PARITY] I don't understand your parity projection %s \n" , 
	    argv[ PARITY_PROJ ] ) ;
    goto memfree ;
  }

  // are we flipping the time direction?
  GLU_bool time_flip = GLU_FALSE ;
  if( are_equal( argv[ TIME_FLIP ] , "true" ) ) {
    time_flip = GLU_TRUE ;
  }

  // get the dimensions
  int mu ;
  char *tok = strtok( (char*)argv[ DIMENSIONS ] , "," ) ;
  Latt.dims[ 0 ] = (int)atoi( tok ) ;
  if( Latt.dims[ 0 ] < 1 ) {
    printf( "[INPUTS] non-sensical lattice dimension %d %d \n" , 
	    0 , Latt.dims[0] ) ;
    goto memfree ;
  }
  for( mu = 1 ; mu < ND ; mu++ ) {
    char *ptok = strtok( NULL , "," ) ;
    if( ptok == NULL ) break ;
    Latt.dims[ mu ] = (int)atoi( ptok ) ;
    if( Latt.dims[mu] < 1 ) {
      printf( "[INPUTS] non-sensical lattice dimension %d %d \n" , 
	      mu , Latt.dims[mu] ) ;
      goto memfree ;
    }
  }

  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;
  if( make_gammas( GAMMAS , basis ) == FAILURE ) {
    goto memfree ;
  }

  // is defined in reader.c
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;

  // allocate projected corr
  proj_corr = allocate_momcorrs( B_CHANNELS , B_CHANNELS , NMOM[0] ) ;

  if( corr == NULL || proj_corr == NULL ) goto memfree ;

  // do the projection
  size_t GSGK ;
#pragma omp parallel for private( GSGK )
  for( GSGK = 0 ; GSGK < ( B_CHANNELS * B_CHANNELS ) ; GSGK++ ) {

    const size_t GSRC = GSGK / B_CHANNELS ;
    const size_t GSNK = GSGK % B_CHANNELS ;

    // loop momenta
    size_t p ;
    for( p = 0 ; p < NMOM[0] ; p++ ) {

      // projections happen here
      const double complex *C = baryon_project( (const struct mcorr**)corr , 
						GAMMAS , momentum ,
						GSRC , GSNK , p ,
						parity_proj ,
						spin_proj ) ;
      
      // poke into proj_corr
      size_t t ;
      if( time_flip == GLU_TRUE ) {
	for( t = 0 ; t < LT ; t++ ) {
	  proj_corr[ GSRC ][ GSNK ].mom[ p ].C[ LT-t-1 ] = C[ t ] ;
	}
      } else {
	for( t = 0 ; t < LT ; t++ ) {
	  proj_corr[ GSRC ][ GSNK ].mom[ p ].C[ t ] = C[ t ] ;
	}
      }
    
      free( (void*)C ) ;
    }
  }

  // write out the correlator
  write_momcorr( argv[ OUTFILE ] , (const struct mcorr **)proj_corr , 
		 momentum , B_CHANNELS , B_CHANNELS , (const int*)NMOM ) ;

 memfree :

  // free the gamma matrices
  free( GAMMAS ) ;

  if( NMOM[0] != 0 ) {
    // free the memory of the read-in correlator
    free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

    // free the memory of the projected correlator
    free_momcorrs( proj_corr , B_CHANNELS , B_CHANNELS , NMOM[0] ) ;
  }

  // free the momentum list
  free( momentum ) ;

  // close the infile
  fclose( infile ) ;

  return 0 ;
}
