/**
   @file correlators.c
   @brief now holds the correlator definitions and IO
 */

#include "common.h"

// allocate the correlator matrix
void
allocate_corrs( struct correlator **corr )
{
  int s ;
  for( s = 0 ; s < NS * NS ; s++ ) {
    corr[ s ] = ( struct correlator* )calloc( NS * NS , sizeof( struct correlator ) ) ;
    int t ;
    for( t = 0 ; t < NS * NS ; t++ ) {
      corr[s][t].C = calloc( L0 , sizeof( double complex ) ) ;
    }
  }
  return ;
}

// free the correlator matrix
void
free_corrs( struct correlator **corr )
{
  int s ;
  for( s = 0 ; s < NS * NS ; s++ ) {
    int t ;
    for( t = 0 ; t < NS * NS ; t++ ) {
      free( corr[ s ][ t ].C ) ;
    }
    free( corr[ s ] ) ;
  }
  free( corr ) ;
  return ;
}

// little convenience funtion
static void
print_convenience( const struct correlator **corr ,
		   const int GSRC ,
		   const int GSNK ) 
{
  int t ;
  for( t = 0 ; t < L0 ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GSRC][GSNK].C[t] ) , cimag( corr[GSRC][GSNK].C[t] ) ) ;
  }
  return ;
}

// debug printing
void
debug_mesons( const char *message , 
	      const struct correlator **corr )
{
  printf( "%s PION\n" , message ) ;
  print_convenience( corr , GAMMA_5 , GAMMA_5 ) ;

  printf( "%s 00\n" , message ) ;
  print_convenience( corr , GAMMA_0 , GAMMA_0 ) ;

  printf( "%s 11\n" , message ) ;
  print_convenience( corr , GAMMA_1 , GAMMA_1 ) ;

  printf( "%s 22\n" , message ) ;
  print_convenience( corr , GAMMA_2 , GAMMA_2 ) ;

  printf( "%s 33\n" , message ) ;
  print_convenience( corr , GAMMA_3 , GAMMA_3 ) ;

  printf( "%s 1010\n" , message ) ;
  print_convenience( corr , 10 , 10 ) ;
  return ;
}
