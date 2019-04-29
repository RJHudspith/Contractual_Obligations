/**
   @file vpffile.c
   @brief VPF momentum-file reader
 */
#include "common.h"

#include "crc32c.h"       // do the crc of the binary data
#include "correlators.h"  // allocate and free corrs
#include "geometry.h"     // init_geom()
#include "GLU_bswap.h"    // byte swaps if necessary

static GLU_bool must_swap = GLU_FALSE ;

struct latt_info Latt ;

// read the binary data
static int
read_data( double *PIdata ,
	   FILE *infile ,
	   const int NMOM )
{
  if( fread( PIdata , sizeof( double ) , NMOM , infile ) != (size_t)NMOM ) {
    fprintf( stderr , "[IO] Unexpected EOF\n" ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_64( NMOM , PIdata ) ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  size_t p ;
  for( p = 0 ; p < (size_t)NMOM ; p++ ) {
    DML_checksum_accum_crc32c( &cksuma , &cksumb , p , 
			       PIdata + p , sizeof( double ) ) ;
  }
  uint32_t cksum[ 2 ] ;
  if( fread( cksum , sizeof( uint32_t ) , 2 , infile ) != 2 ) {
    fprintf( stderr , "[IO] Cannot find checksums\n" ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_32( 2 , cksum ) ;
  if( cksuma != cksum[0] || cksumb != cksum[1] ) {
    fprintf( stderr , "[IO] Checksum mismatch!! \n" ) ;
    fprintf( stderr , "[IO] Computed %x %x :: Read %x %x \n" ,
	     cksuma , cksumb , cksum[0] , cksum[1] ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// quick little accessor
static int
FREAD32( uint32_t *data , const size_t size , FILE *infile ) {
  if( fread( data , sizeof( uint32_t ) , size , infile ) != size ) {
    fprintf( stderr , "[IO] FREAD32 failure \n" ) ;
    return FAILURE ;
  }
  if( must_swap ) bswap_32( size , data ) ;
  return SUCCESS ;
}

// write it out
static void
write_PIdata( const int **momentum ,
	      const double *PIdata ,
	      const size_t NMOM ,
	      const size_t DIMS )
{
  size_t p ;
  fprintf( stdout , "\n[MOMS] outputting available %zu-momenta ... \n\n" , 
	   DIMS ) ;
  for( p = 0 ; p < NMOM ; p++ ) {
    register double p2 = 0.0 ;
    size_t mu ;
    fprintf( stdout , "[MOMS] %zu :: (" , p ) ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      const double cache = 2.0 * sin( 0.5 * momentum[ p ][ mu ] * TWOPI / (double)Latt.dims[ mu ] ) ;
      p2 += cache * cache ;
      fprintf( stdout , " %d " , momentum[ p ][ mu ] ) ;
    }
    fprintf( stdout , ") %f %f \n" , p2 , PIdata[ p ] ) ;
  }
  fprintf( stdout , "\n" ) ;
  return ;
}

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc != 3 ) {
    return fprintf( stdout , "[VPFREAD] "
		    "usage ./VPF {vpf_file} lx,ly,lz,lt \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[IO] File %s does not exist\n" , argv[1] ) ;
    return FAILURE ;
  }

  // set up lattice geometry
  {
    char *tok1 = strtok( (char*)argv[2] , "," ) ;
    Latt.dims[0] = (size_t)atoi( tok1 ) ;

    int mu = 1 ;
    while( ( tok1 = strtok( NULL , "," ) ) != NULL ) {
      Latt.dims[ mu ] = (size_t)atoi( tok1 ) ;
      mu++ ;
    }
    if( mu < ND || mu > ND ) {
      fprintf( stderr , "[IO] Incorrect lattice dimensions ::"
	       " Compiled %d, Given %d, \n" , ND , mu ) ;
      return FAILURE ;
    }
  }
  init_geom( ) ;

  // read the file
  uint32_t magic[1] , NMOM[1] = { 0 } ;

  // PIdata and 
  double *PIdata = NULL ;
  int **momentum = NULL ;

  GLU_bool failure = GLU_FALSE ;

  if( fread( magic , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    goto memfree ;
  }

  // check the magic number, tells us the edianness
  if( magic[0] != VPF_MAGIC ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != VPF_MAGIC) {
      fprintf( stderr , "[IO] Magic number read failure\n" ) ;
      goto memfree ;
    }
    must_swap = GLU_TRUE ;
  }

  // read the length of the momentum list
  if( FREAD32( NMOM , 1 , infile ) == FAILURE ) goto memfree ;

  // this is a logic error
  if( NMOM[0] < 1 ) goto memfree ;

  momentum = malloc( NMOM[0] * sizeof( int* ) ) ;

  int p ;
  for( p = 0 ; p < (int)NMOM[0] ; p++ ) {
    momentum[ p ] = malloc( ( ND ) * sizeof( int ) ) ;
    uint32_t n[ ND + 1 ] ;
    if( FREAD32( n , ND + 1 , infile ) == FAILURE ) {
      failure = GLU_TRUE ;
    }
    if( n[ 0 ] != ND ) {
      fprintf( stderr , "[MOMLIST] %d should be %d \n" , n[ 0 ] , ND ) ;
      failure = GLU_TRUE ;
    }
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      momentum[ p ][ mu ] = (int)n[ 1 + mu ] ;
    }
  }
  if( failure == GLU_TRUE ) {
    goto memfree ;
  }

  // read in the momentum list size again 
  uint32_t TNMOM[ 1 ] ;
  if( FREAD32( TNMOM , 1 , infile ) == FAILURE ) goto memfree ;
  if( TNMOM[ 0 ] != NMOM[ 0 ] ) {
    fprintf( stderr , "[MOMLIST] length mismatch %d %d \n" , 
	     NMOM[0] , TNMOM[0] ) ;
    goto memfree ;
  }

  // read the data
  PIdata = malloc( NMOM[0] * sizeof( double ) ) ;
  if( read_data( PIdata , infile , NMOM[0] ) == FAILURE ) goto memfree ;

  write_PIdata( (const int**)momentum , PIdata , NMOM[0] , ND ) ;

 memfree :

  free( PIdata ) ;

  // free the momentum list
  if( momentum != NULL ) {
    for( p = 0 ; p < (int)NMOM[ 0 ] ; p++ ) {
      free( momentum[ p ] ) ;
    }
  }
  free( momentum ) ;

  fclose( infile ) ;

  return 0 ;
}

