
#include "common.h"

#include "crc32.h"        // do the crc of the binary data
#include "correlators.h"  // allocate and free corrs
#include "geometry.h"     // init_geom()
#include "GLU_bswap.h"    // byte swaps if necessary

static GLU_bool must_swap = GLU_FALSE ;

struct latt_info Latt ;

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 2 ) {
    return printf( "usage ./MESONS {correlator file} GSRC,GSNK ... \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return -1 ;
  }

  uint32_t magic[1] , NGSRC[1] , NGSNK[1] ;

  if( fread( magic , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    return FAILURE ;
  }

  // check the magic number, tells us the edianness
  if( magic[0] != 67798233 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67798233 ) {
      printf( "Magic number read failure\n" ) ;
      return FAILURE ;
    }
    must_swap = GLU_TRUE ;
  }

  if( fread( NGSRC , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , NGSRC ) ;

  if( fread( NGSNK , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , NGSNK ) ;

  uint32_t LT[1] ;
  if( fread( LT , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 1 , LT ) ;
  Latt.dims[0] = 0 ;
  Latt.dims[1] = 0 ;
  Latt.dims[2] = 0 ;
  Latt.dims[ ND - 1 ] = (int)LT[ 0 ] ;

  init_geom( ) ;

  struct correlator **corr  = malloc( NS * NS * sizeof( struct correlator* ) ) ;
  allocate_corrs( corr ) ;

  // read the correlator
  int GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < NGSRC[0] ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < NGSNK[0] ; GSNK++ ) {
      if( GSRC == 0 && GSNK == 0 ) {
	// read the correlator
	if( fread( corr[ GSRC ][ GSNK ].C , sizeof( double complex ) , 
		   L0 , infile ) != L0 ) {
	  return FAILURE ;
	}
	if( must_swap ) bswap_64( L0 * 2 , corr[ GSRC ][ GSNK ].C ) ;
      } else {
	uint32_t LT[1] ;
	if( fread( LT , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
	  return FAILURE ;
	}
	if( must_swap ) bswap_32( 1 , LT ) ;
	if( fread( corr[ GSRC ][ GSNK ].C , sizeof( double complex ) , 
		   L0 , infile ) != L0 ) {
	  return FAILURE ;
	}
	if( must_swap ) bswap_64( L0 * 2 , corr[ GSRC ][ GSNK ].C ) ;
      }
      //
    }
  }

  // loop the ones we want
  int i ;
  for( i = 2 ; i < ( argc ) ; i++ ) {
    // tokenize argv into the correlators people want
    char *tok1 = strtok( (char*)argv[i] , "," ) ;
    char *tok2 = strtok( NULL , "," ) ;

    if( tok1 == NULL || tok2 == NULL ) break ;

    const int idx1 = atoi( tok1 ) ;
    const int idx2 = atoi( tok2 ) ;

    printf( "Correlator [ %d %d ] \n" , idx1 , idx2 ) ;
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      printf( "%d %1.12e %1.12e\n" , t ,
	      creal( corr[ idx1 ][ idx2 ].C[ t ] ) ,
	      cimag( corr[ idx1 ][ idx2 ].C[ t ] ) ) ;
    }
    //
  }
 
  free_corrs( corr ) ;

  fclose( infile ) ;

  return 0 ;
}
