/**
   @file correlators.c
   @brief now holds the correlator definitions and IO
 */

#include "common.h"

#include "crc32.h"        // crc 
#include "cut_output.h"   // for write_mom_veclist

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

// allocate the correlator matrix
struct correlator **
allocate_corrs( const int NSRC , 
		const int NSNK )
{
  struct correlator **corr = (struct correlator**)malloc( NSRC * sizeof( struct correlator* ) ) ; 
  int GSRC ;
  for( GSRC = 0 ; GSRC < NSRC ; GSRC++ ) {
    corr[ GSRC ] = ( struct correlator* )malloc( NSNK * sizeof( struct correlator ) ) ;
    int GSNK ;
    for( GSNK = 0 ; GSNK < NSNK ; GSNK++ ) {
      corr[ GSRC ][ GSNK ].C = calloc( L0 , sizeof( double complex ) ) ;
    }
  }
  return corr ;
}

// allocation of dispersion relation correlation function
struct mcorr **
allocate_momcorrs( const int length1 , 
		   const int length2 ,
		   const int nmom )
{
  struct mcorr **mcorr = malloc( length1 * sizeof( struct mcorr* ) ) ;
  int i ;
  for( i = 0 ; i < length1 ; i++ ) {
    mcorr[ i ] = malloc( length2 * sizeof( struct mcorr ) ) ;
    int j ;
    for( j = 0 ; j < length2 ; j++ ) {
      mcorr[ i ][ j ].mom = malloc( nmom * sizeof( struct correlator ) ) ;
      int p ;
      for( p = 0 ; p < nmom ; p++ ) {
	mcorr[ i ][ j ].mom[ p ].C = malloc( L0 * sizeof( double complex ) ) ;
      }
    }
  }
  return mcorr ;
}

// free the correlator matrix
void
free_corrs( struct correlator **corr , 
	    const int NSRC ,
	    const int NSNK )
{
  int s ;
  for( s = 0 ; s < NSRC ; s++ ) {
    int t ;
    for( t = 0 ; t < NSNK ; t++ ) {
      free( corr[ s ][ t ].C ) ;
    }
    free( corr[ s ] ) ;
  }
  free( corr ) ;
  return ;
}

// momcorr freer
void
free_momcorrs( struct mcorr **mcorr , 
	       const int length1 ,
	       const int length2 ,
	       const int nmom ) 
{
  int i ;
  for( i = 0 ; i < length1 ; i++ ) {
    int j ;
    for( j = 0 ; j < length2 ; j++ ) {
      int p ;
      for( p = 0 ; p < nmom ; p++ ) {
	free( mcorr[ i ][ j ].mom[ p ].C ) ;
      }
      free( mcorr[ i ][ j ].mom ) ;
    }
    free( mcorr[ i ] ) ;
  }
  free( mcorr ) ;
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

// debug printing for baryons
void
debug_baryons( const char *message , 
	      const struct correlator **corr )
{
  printf( "%s OCTETT\n" , message ) ;
  print_convenience( corr , 5 , 0 ) ;

  printf( "%s DECUPLETT G0\n" , message ) ;
  print_convenience( corr , 0 , 0 ) ;

  printf( "%s DECUPLETT G1\n" , message ) ;
  print_convenience( corr , 1 , 0 ) ;

  printf( "%s DECUPLETT G2\n" , message ) ;
  print_convenience( corr , 2 , 0 ) ;

  printf( "%s DECUPLETT G3\n" , message ) ;
  print_convenience( corr , 3 , 0 ) ;

  printf( "%s DECUPLETT G4\n" , message ) ;
  print_convenience( corr , 4 , 0 ) ;

  return ;
}

// write the full correlator matrix
void
write_correlators( const char *outfile ,
		   const struct correlator **corr ,
		   const int NSRC ,
		   const int NSNK )
{
  printf( "[IO] writing correlation matrix to %s \n" , outfile ) ;

  FILE *output_file = fopen( outfile , "wb" ) ;

  uint32_t magic[ 1 ] = { 67798233 } ; // THIS SPELLS COR! in ascii

  fwrite( magic , sizeof( uint32_t ) , 1 , output_file ) ;

  uint32_t NGSRC[ 1 ] = { NSRC } ;
  uint32_t NGSNK[ 1 ] = { NSNK } ;

  fwrite( NGSRC , sizeof( uint32_t ) , 1 , output_file ) ;
  fwrite( NGSNK , sizeof( uint32_t ) , 1 , output_file ) ;

  uint32_t LT[ 1 ] = { L0 } , cksuma = 0 , cksumb = 0 ;
  int GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < NSRC ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < NSNK ; GSNK++ ) {
      fwrite( LT , sizeof( uint32_t ) , 1 , output_file ) ;
      fwrite( corr[GSRC][GSNK].C , sizeof( double complex ) , L0 , output_file ) ; 
      // accumulate a checksum
      DML_checksum_accum( &cksuma , &cksumb , GSNK + NSNK * GSRC , 
			  (char*)corr[GSRC][GSNK].C , 
			  sizeof( double complex ) * L0 ) ;
    }
  }
  // write out both checksums
  uint32_t csum[ 2 ] = { cksuma , cksumb } ;
  fwrite( csum , sizeof( uint32_t ) , 2 , output_file ) ;

  fclose( output_file ) ;

  return ;
}

// write the full correlator matrix
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const int NSRC ,
	       const int NSNK ,
	       const int *nmom )
{
  printf( "[IO] writing correlation matrix to %s \n" , outfile ) ;

  FILE *output_file = fopen( outfile , "wb" ) ;

  uint32_t magic[ 1 ] = { 67798233 } ; // THIS SPELLS COR! in ascii

  fwrite( magic , sizeof( uint32_t ) , 1 , output_file ) ;

  write_mom_veclist( output_file , nmom , list , ND-1 ) ;

  uint32_t NMOM[ 1 ] = { nmom[0] } ;
  
  fwrite( NMOM , sizeof( uint32_t ) , 1 , output_file ) ;

  uint32_t LT[ 1 ] = { L0 } , cksuma = 0 , cksumb = 0 ;

  int p ;
  for( p = 0 ; p < nmom[0] ; p++ ) {
    
    uint32_t NGSRC[ 1 ] = { NSRC } ;
    uint32_t NGSNK[ 1 ] = { NSNK } ;
    
    fwrite( NGSRC , sizeof( uint32_t ) , 1 , output_file ) ;
    fwrite( NGSNK , sizeof( uint32_t ) , 1 , output_file ) ;
    
    int GSRC , GSNK ;
    for( GSRC = 0 ; GSRC < NSRC ; GSRC++ ) {
      for( GSNK = 0 ; GSNK < NSNK ; GSNK++ ) {
	fwrite( LT , sizeof( uint32_t ) , 1 , output_file ) ;
	fwrite( corr[GSRC][GSNK].mom[p].C , sizeof( double complex ) , L0 , output_file ) ; 
	// accumulate a checksum
	DML_checksum_accum( &cksuma , &cksumb , GSNK + NSNK * GSRC , 
			    (char*)corr[GSRC][GSNK].mom[p].C , 
			    sizeof( double complex ) * L0 ) ;
      }
    }
  }
  // write out both checksums
  uint32_t csum[ 2 ] = { cksuma , cksumb } ;
  fwrite( csum , sizeof( uint32_t ) , 2 , output_file ) ;

  fclose( output_file ) ;

  return ;
}
