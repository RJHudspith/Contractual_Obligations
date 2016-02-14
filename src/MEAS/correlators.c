/**
   @file correlators.c
   @brief now holds the correlator definitions and IO
 */
#include "common.h"

#include "crc32c.h"       // crc 
#include "correlators.h"  // so that I can alphabetise
#include "cut_output.h"   // for write_mom_veclist

// little convenience funtion
static void
print_convenience( const struct mcorr **corr ,
		   const size_t GSRC ,
		   const size_t GSNK ) 
{
  int t ;
  for( t = 0 ; t < LT ; t++ ) {
    printf( "%d %e %e \n" , t , creal( corr[GSRC][GSNK].mom[0].C[t] ) , 
	    cimag( corr[GSRC][GSNK].mom[0].C[t] ) ) ;
  }
  return ;
}

// allocation of dispersion relation correlation function
struct mcorr **
allocate_momcorrs( const size_t length1 , 
		   const size_t length2 ,
		   const size_t nmom )
{
  // check if stuff is non-zero
  if( length1 == 0 || length2 == 0 || nmom == 0 ) {
    printf( "[IO] corr lengths are 0 :: ( %zu , %zu , %zu ) \n" ,
	    length1 , length2 , nmom ) ;
    return NULL ;
  }
  struct mcorr **mcorr = malloc( length1 * sizeof( struct mcorr* ) ) ;
  size_t i , j , p ;
  for( i = 0 ; i < length1 ; i++ ) {
    mcorr[ i ] = malloc( length2 * sizeof( struct mcorr ) ) ;
    for( j = 0 ; j < length2 ; j++ ) {
      mcorr[ i ][ j ].mom = malloc( nmom * sizeof( struct correlator ) ) ;
      for( p = 0 ; p < nmom ; p++ ) {
	mcorr[ i ][ j ].mom[ p ].C = malloc( LT * sizeof( double complex ) ) ;
      }
    }
  }
  return mcorr ;
}

// momcorr freer
void
free_momcorrs( struct mcorr **mcorr , 
	       const size_t length1 ,
	       const size_t length2 ,
	       const size_t nmom ) 
{
  size_t i , j , p ;
  for( i = 0 ; i < length1 ; i++ ) {
    for( j = 0 ; j < length2 ; j++ ) {
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
	      const struct mcorr **corr )
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
	       const struct mcorr **corr )
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

// baryon writers
void
write_baryon( struct mcorr **corr , 
	      const struct veclist *list ,
	      const int NMOM[ 1 ] ,
	      const GLU_bool is_wall ,
	      const char *outfile ,
	      const char *type )
{
  // write out the correlator
  char outstr[ 256 ] , wall[ 8 ] = "" ;
  if( is_wall == GLU_TRUE ) {
    sprintf( wall , ".WW" ) ;
  } 

  // write out whatever type we are using
  sprintf( outstr , "%s.%s%s" , outfile , type , wall ) ;
  write_momcorr( outstr , (const struct mcorr**)corr , 
		 list , B_CHANNELS * B_CHANNELS , NSNS , NMOM ) ;

  return ;
}

// write the full correlator matrix
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const size_t NSRC ,
	       const size_t NSNK ,
	       const int *nmom )
{
  printf( "[IO] writing correlation matrix to %s \n" , outfile ) ;

  FILE *output_file = fopen( outfile , "wb" ) ;

  uint32_t magic[ 1 ] = { 67798233 } ; // THIS SPELLS COR! in ascii

  fwrite( magic , sizeof( uint32_t ) , 1 , output_file ) ;

  write_mom_veclist( output_file , nmom , list , ND-1 ) ;

  uint32_t NMOM[ 1 ] = { nmom[0] } ;
  
  fwrite( NMOM , sizeof( uint32_t ) , 1 , output_file ) ;

  uint32_t L0[ 1 ] = { LT } , cksuma = 0 , cksumb = 0 ;

  size_t p ;
  for( p = 0 ; p < nmom[0] ; p++ ) {
    
    uint32_t NGSRC[ 1 ] = { (uint32_t)NSRC } ;
    uint32_t NGSNK[ 1 ] = { (uint32_t)NSNK } ;
    
    fwrite( NGSRC , sizeof( uint32_t ) , 1 , output_file ) ;
    fwrite( NGSNK , sizeof( uint32_t ) , 1 , output_file ) ;

    size_t GSRC , GSNK ;
    for( GSRC = 0 ; GSRC < NSRC ; GSRC++ ) {
      for( GSNK = 0 ; GSNK < NSNK ; GSNK++ ) {
	fwrite( L0 , sizeof( uint32_t ) , 1 , output_file ) ;
	fwrite( corr[GSRC][GSNK].mom[p].C , sizeof( double complex ) , LT , output_file ) ; 
	// accumulate the newer, fancier checksum
	DML_checksum_accum_crc32c( &cksuma , &cksumb ,
				   GSNK + NSNK * GSRC , 
				   corr[GSRC][GSNK].mom[p].C , 
				   sizeof( double complex ) * LT ) ;
      }
    }
  }
  // write out both checksums
  uint32_t csum[ 2 ] = { cksuma , cksumb } ;
  fwrite( csum , sizeof( uint32_t ) , 2 , output_file ) ;

  fclose( output_file ) ;

  return ;
}
