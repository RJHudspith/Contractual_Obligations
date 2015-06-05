/**
   @file mesonfile.c
   @brief little correlation file reader
 */
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
    return FAILURE ;
  }
  if( must_swap ) bswap_64( L0 * 2 , corr ) ;

  // accumulate the checksums
  DML_checksum_accum( cksuma , cksumb , rank , 
		      (char*)corr, 
		      sizeof( double complex ) * L0 ) ;
  return SUCCESS ;
}

// quick little accessor
static int
FREAD32( uint32_t *data , const int size , FILE *infile ) {
  if( fread( data , sizeof( uint32_t ) , size , infile ) != size ) {
    printf( "[IO] FREAD32 failure \n" ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_32( size , data ) ;
  return SUCCESS ;
}

// read the momentum-correlator
static int
read_momcorr( struct mcorr **corr ,
	      FILE *infile ,
	      uint32_t *cksuma ,
	      uint32_t *cksumb ,
	      const uint32_t NGSRC[ 1 ] ,
	      const uint32_t NGSNK[ 1 ] ,
	      const uint32_t LT[ 1 ] ,
	      const int p )
{
  if( p != 0 ) {
    uint32_t tNGSRC[ 1 ] , tNGSNK[ 1 ] , tLT[ 1 ] ;
    if( FREAD32( tNGSRC , 1 , infile ) == FAILURE ) return FAILURE ;
    if( FREAD32( tNGSNK , 1 , infile ) == FAILURE ) return FAILURE ;
    if( FREAD32( tLT , 1 , infile ) == FAILURE ) return FAILURE ;
    if( tNGSRC[ 0 ] != NGSRC[ 0 ] || tNGSNK[ 0 ] != NGSNK[ 0 ] || tLT[ 0 ] != LT[ 0 ] ) {
      printf( "[IO] NGSRC , NGSNK , LT misread \n" ) ;
      printf( "[IO] %d %d %d \n" , tNGSRC[ 0 ] , tNGSNK[ 0 ] , tLT[ 0 ] ) ;
      return FAILURE ;
    }
  }

  int GSRC , GSNK ;
  for( GSRC = 0 ; GSRC < (int)NGSRC[0] ; GSRC++ ) {
    for( GSNK = 0 ; GSNK < (int)NGSNK[0] ; GSNK++ ) {
      if( ( GSRC == 0 ) && ( GSNK == 0 ) ) {
      } else {
	// read the timeslice stuff
	uint32_t LT[1] ;
	if( FREAD32( LT , 1 , infile ) == FAILURE ) return FAILURE ;
	if( (int)LT[0] != L0 ) { 
	  printf( "[IO] LT Read failure %d %d \n" , (int)LT[0] , L0 ) ; 
	  return FAILURE ; 
	}
      }
      if( read_corr( corr[ GSRC ][ GSNK ].mom[ p ].C , cksuma , cksumb , infile ,
		     GSNK + NGSNK[0] * GSRC ) == FAILURE ) {
	printf( "[IO] corr Read failure \n" ) ;
	return FAILURE ;
      }
      //
    }
  }
  return SUCCESS ;
}

// finds the desired mom
static int
find_desired_mom( const int **momentum , 
		  const int *moms , 
		  const int NMOM )
{
  int i ;
  for( i = 0 ; i < NMOM ; i++ ) {
    int mu , matches = 0 ;
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      if( momentum[ i ][ mu ] == moms[ mu ] ) {
	matches++ ;
      }
    }
    if( matches == ND-1 ) return i ;
  }
  return FAILURE ;
}

// this is for you, Anthony
static void
write_momlist( const int **momentum ,
	       const int NMOM )
{
  int p ;
  printf( "\n[MOMS] outputting available %d-momenta ... \n\n" , ND - 1 ) ;
  for( p = 0 ; p < NMOM ; p++ ) {
    int mu ;
    printf( "[MOMS] %d :: (" , p ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      printf( " %d " , momentum[ p ][ mu ] ) ;
    }
    printf( ") \n" ) ;
  }
  printf( "\n" ) ;
  return ;
}

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 2 ) {
    return printf( "usage ./MESONS {correlator file} GSRC,GSNK,px,py,pz ... \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return -1 ;
  }

  uint32_t magic[1] , NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL ;
  int **momentum = NULL ;

  // checksums
  uint32_t cksuma = 0 , cksumb = 0 , csum[ 2 ] = { 0 , 0 } ; 

  // number of correlators printed to stdout
  int corrs_written = 0 ;

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

  // read the length of the momentum list
  if( FREAD32( NMOM , 1 , infile ) == FAILURE ) return FAILURE ;

  momentum = malloc( NMOM[0] * sizeof( int* ) ) ;

  GLU_bool failure = GLU_FALSE ;
  int p ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    uint32_t n[ ND ] ;
    if( FREAD32( n , ND , infile ) == FAILURE ) failure = GLU_TRUE ;
    momentum[ p ] = malloc( ( ND - 1 ) * sizeof( int ) ) ;
    if( n[ 0 ] != ND-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , ND-1 ) ;
      failure = GLU_TRUE ;
    }
    int mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      momentum[ p ][ mu ] = (int)n[ 1 + mu ] ;
    }
  }
  if( failure == GLU_TRUE ) goto memfree ;

  // read in the momentum list size again 
  uint32_t TNMOM[ 1 ] ;
  if( FREAD32( TNMOM , 1 , infile ) == FAILURE ) goto memfree ;
  if( TNMOM[ 0 ] != NMOM[ 0 ] ) {
    printf( "[MOMLIST] length mismatch %d %d \n" , NMOM[0] , TNMOM[0] ) ;
    goto memfree ;
  }

  if( FREAD32( NGSRC , 1 , infile ) == FAILURE ) goto memfree ;
  if( FREAD32( NGSNK , 1 , infile ) == FAILURE ) goto memfree ;

  // read in an LT
  uint32_t LT[ 1 ] ;
  if( FREAD32( LT , 1 , infile ) == FAILURE ) goto memfree ;

  Latt.dims[0] = 0 ;
  Latt.dims[1] = 0 ;
  Latt.dims[2] = 0 ;
  Latt.dims[ ND - 1 ] = (int)LT[ 0 ] ;

  init_geom( ) ;

  // allocate the correlators after we have got the momentum information
  corr = allocate_momcorrs( NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // read the correlator
  for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
    if( read_momcorr( corr , infile , &cksuma , &cksumb ,
		      NGSRC , NGSNK , LT , p ) == FAILURE ) {
      goto memfree ;
    }
  }

  printf( "[IO] All correlators read \n" ) ;

  // check our checksums
  if( FREAD32( csum , 2 , infile ) == FAILURE ) goto memfree ;
  if( csum[0] != cksuma || csum[1] != cksumb ) {
    printf( "Mismatched checksums ! %x %x %x %x\n" , csum[0] , csum[1] , cksuma , cksumb ) ;
    goto memfree ;
  } 

  printf( "[CHECKSUM] both checksums passed \n\n" ) ;

  // loop the ones we want
  int i ;
  for( i = 2 ; i < ( argc ) ; i++ ) {
    // tokenize argv into the correlators people want
    char *tok1 = strtok( (char*)argv[i] , "," ) ;
    if( tok1 == NULL ) break ;
    const int idx1 = atoi( tok1 ) ;
    if( idx1 >= NGSRC[0] || idx1 < 0 ) { 
      printf( "[Momcorr] Non-sensical source index %d \n" , idx1 ) ;
      break ;
    } 
    char *tok2 = strtok( NULL , "," ) ;
    if( tok2 == NULL ) break ;
    const int idx2 = atoi( tok2 ) ;
    if( idx2 >= NGSNK[0] || idx2 < 0 ) { 
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
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      printf( "CORR %d %1.12e %1.12e\n" , t ,
	      creal( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ,
	      cimag( corr[ idx1 ][ idx2 ].mom[ matchmom ].C[ t ] ) ) ;
    }
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

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
    free( momentum[ p ] ) ;
  }
  free( momentum ) ;

  fclose( infile ) ;

  return 0 ;
}
