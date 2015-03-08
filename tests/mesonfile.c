
#include "common.h"

#include "crc32.h"        // do the crc of the binary data
#include "correlators.h"  // allocate and free corrs
#include "geometry.h"     // init_geom()
#include "GLU_bswap.h"    // byte swaps if necessary

static GLU_bool must_swap = GLU_FALSE ;

struct latt_info Latt ;

// read the correlator for a particular source/sink gamma
static int
read_corr( double complex *corr ,
	   uint32_t *cksuma , 
	   uint32_t *cksumb ,
	   FILE *infile ,
	   const int rank )
{
  if( fread( corr , sizeof( double complex ) , L0 , infile ) != L0 ) {
    printf( "corr reading?\n" ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_64( L0 * 2 , corr ) ;

  // accumulate the checksums
  DML_checksum_accum( cksuma , cksumb , rank , 
		      (char*)corr, 
		      sizeof( double complex ) * L0 ) ;
  return SUCCESS ;
}

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

  // read in an LT
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
  uint32_t cksuma = 0 , cksumb = 0 ; 
  int GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < (int)NGSRC[0] ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < (int)NGSNK[0] ; GSNK++ ) {
      if( ( GSRC == 0 ) && ( GSNK == 0 ) ) {
      } else {
	// read the timeslice stuff
	uint32_t LT[1] ;
	if( fread( LT , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
	  return FAILURE ;
	}
	if( must_swap ) bswap_32( 1 , LT ) ;
	if( (int)LT[0] != L0 ) { 
	  printf( "LT Read failure %d %d \n" , (int)LT[0] , L0 ) ; 
	  return FAILURE ; 
	}
      }
      if( read_corr( corr[ GSRC ][ GSNK ].C , &cksuma , &cksumb , infile ,
		     GSNK + NSNS * GSRC ) == FAILURE ) {
	printf( "Read failure \n" ) ;
	return FAILURE ;
      }
      //
    }
  }

  // check our checksums
  uint32_t csum[ 2 ] ; 
  if( fread( csum , sizeof( uint32_t ) , 2 , infile ) != 1 ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 2 , csum ) ;
  if( csum[0] != cksuma || csum[1] != cksumb ) {
    printf( "Mismatched checksums ! %x %x %x %x\n" , csum[0] , csum[1] , cksuma , cksumb ) ;
    return FAILURE ;
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