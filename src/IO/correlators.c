/**
   @file correlators.c
   @brief now holds the correlator definitions and IO
 */
#include "common.h"

#include "crc32c.h"       // crc 
#include "correlators.h"  // so that I can alphabetise
#include "cut_output.h"   // for write_mom_veclist

// allocation of dispersion relation correlation function
struct mcorr **
allocate_momcorrs( const size_t length1 , 
		   const size_t length2 ,
		   const size_t nmom )
{
  // check if stuff is non-zero
  if( length1 == 0 || length2 == 0 || nmom == 0 ) {
    fprintf( stderr , "[IO] corr lengths are 0 :: ( %zu , %zu , %zu ) \n" ,
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

// write the full correlator matrix
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const size_t NSRC ,
	       const size_t NSNK ,
	       const int *nmom , 
	       const char *type )
{
  // write out the correlator
  char outstr[ 256 ] ;
  if( !strcmp( type , "" ) ) { 
    sprintf( outstr , "%s" , outfile ) ;
  } else {
    sprintf( outstr , "%s.%s" , outfile , type ) ;
  }

  fprintf( stdout , "[IO] writing correlation matrix to %s \n" , outstr ) ;

  FILE *output_file = fopen( outstr , "wb" ) ;

  uint32_t magic[ 1 ] = { CORR_MAGIC } ; // THIS SPELLS CORR in ascii

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
				   p + nmom[0] * ( GSNK + NSNK * GSRC ) , 
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
