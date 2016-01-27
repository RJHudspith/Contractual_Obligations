/**
   @file reader.c
   @brief correlator file reader
 */
#include "common.h"

#include "crc32c.h"        // do the crc of the binary data
#include "correlators.h"  // allocate and free corrs
#include "geometry.h"     // init_geom()
#include "GLU_bswap.h"    // byte swaps if necessary

// do we have to byte swap?
static GLU_bool must_swap = GLU_FALSE ;

// read the correlator for a particular source/sink gamma
static int
read_corr( double complex *corr ,
	   uint32_t *cksuma , 
	   uint32_t *cksumb ,
	   FILE *infile ,
	   const int rank )
{
  if( fread( corr , sizeof( double complex ) , LT , infile ) != LT ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_64( LT * 2 , corr ) ;

  // accumulate the checksums
  DML_checksum_accum_crc32c( cksuma , cksumb , rank , 
			     corr , 
			     sizeof( double complex ) * LT ) ;
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
	      const uint32_t L0[ 1 ] ,
	      const int p )
{
  if( p != 0 ) {
    uint32_t tNGSRC[ 1 ] , tNGSNK[ 1 ] , tLT[ 1 ] ;
    if( FREAD32( tNGSRC , 1 , infile ) == FAILURE ) return FAILURE ;
    if( FREAD32( tNGSNK , 1 , infile ) == FAILURE ) return FAILURE ;
    if( FREAD32( tLT , 1 , infile ) == FAILURE ) return FAILURE ;
    if( tNGSRC[ 0 ] != NGSRC[ 0 ] || tNGSNK[ 0 ] != NGSNK[ 0 ] || tLT[ 0 ] != L0[ 0 ] ) {
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
	uint32_t L0[1] ;
	if( FREAD32( L0 , 1 , infile ) == FAILURE ) return FAILURE ;
	if( (int)L0[0] != LT ) { 
	  printf( "[IO] LT Read failure %d %d \n" , (int)L0[0] , LT ) ; 
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
size_t
find_desired_mom( const struct veclist *momentum , 
		  const int *moms , 
		  const int NMOM )
{
  size_t i ;
  for( i = 0 ; i < NMOM ; i++ ) {
    size_t mu , matches = 0 ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      if( momentum[ i ].MOM[ mu ] != moms[ mu ] ) break ;
      matches++ ;
    }
    if( matches == ND-1 ) return i ;
  }
  return FAILURE ;
}

// this is for you, Anthony
void
write_momlist( const struct veclist *momentum ,
	       const int NMOM )
{
  int p ;
  printf( "\n[MOMS] outputting available %d-momenta ... \n\n" , ND - 1 ) ;
  for( p = 0 ; p < NMOM ; p++ ) {
    int mu ;
    printf( "[MOMS] %d :: (" , p ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      printf( " %d " , momentum[ p ].MOM[ mu ]  ) ;
    }
    printf( ") \n" ) ;
  }
  printf( "\n" ) ;
  return ;
}

// allocates and packs mcorr array
struct mcorr**
process_file( struct veclist **momentum ,
	      FILE *infile ,
	      uint32_t NGSRC[1] ,
	      uint32_t NGSNK[1] ,
	      uint32_t NMOM[1] )
{
  // magic number
  uint32_t magic[1] ;

  // checksums
  uint32_t cksuma = 0 , cksumb = 0 , csum[ 2 ] = { 0 , 0 } ; 

  // correlators
  struct mcorr **corr = NULL ;

  // read the magic number
  if( fread( magic , sizeof( uint32_t ) , 1 , infile ) != 1 ) {
    return NULL ;
  }

  // check the magic number, tells us the edianness
  if( magic[0] != 67798233 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67798233 ) {
      printf( "Magic number read failure\n" ) ;
      return NULL ;
    }
    must_swap = GLU_TRUE ;
  }

  // read the length of the momentum list
  if( FREAD32( NMOM , 1 , infile ) == FAILURE ) return NULL ;

  *momentum = malloc( NMOM[0] * sizeof( struct veclist ) ) ;

  GLU_bool failure = GLU_FALSE ;
  size_t p ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    uint32_t n[ ND ] ;
    if( FREAD32( n , ND , infile ) == FAILURE ) failure = GLU_TRUE ;
    if( n[ 0 ] != ND-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , ND-1 ) ;
      failure = GLU_TRUE ;
    }
    size_t mu ;
    (*momentum)[ p ].nsq = 0 ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      (*momentum)[ p ].MOM[ mu ] = (int)n[ 1 + mu ] ;
      (*momentum)[ p ].nsq += n[ 1 + mu ] * n[ 1 + mu ] ;
    }
    (*momentum)[ p ].MOM[ mu ] = 0 ; // set the t-direction to 0
  }
  if( failure == GLU_TRUE ) return NULL ;

  // read in the momentum list size again 
  uint32_t TNMOM[ 1 ] ;
  if( FREAD32( TNMOM , 1 , infile ) == FAILURE ) return NULL ;
  if( TNMOM[ 0 ] != NMOM[ 0 ] ) {
    printf( "[MOMLIST] length mismatch %d %d \n" , NMOM[0] , TNMOM[0] ) ;
    return NULL ;
  }

  if( FREAD32( NGSRC , 1 , infile ) == FAILURE ) return NULL ;
  if( FREAD32( NGSNK , 1 , infile ) == FAILURE ) return NULL ;

  // read in an LT
  uint32_t L0[ 1 ] ;
  if( FREAD32( L0 , 1 , infile ) == FAILURE ) return NULL ;

  Latt.dims[ ND - 1 ] = (int)L0[ 0 ] ;

  init_geom( ) ;

  // allocate the correlators after we have got the momentum information
  corr = allocate_momcorrs( NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // read the correlator
  for( p = 0 ; p < NMOM[ 0 ] ; p++ ) {
    if( read_momcorr( corr , infile , &cksuma , &cksumb ,
		      NGSRC , NGSNK , L0 , p ) == FAILURE ) {
      return NULL ;
    }
  }

  printf( "[IO] All correlators read \n" ) ;

  // check our checksums
  // at the moment this is a soft warning as I changed the type of checksum we 
  // use to one that has an intrinsic
  // TODO :: change it to a full exit in the near future - J
  if( FREAD32( csum , 2 , infile ) == FAILURE ) return NULL ;
  if( csum[0] != cksuma || csum[1] != cksumb ) {
    printf( "[CHECKSUM] Mismatched checksums ! %x %x %x %x\n" , csum[0] , csum[1] , cksuma , cksumb ) ;
  } else {
    printf( "[CHECKSUM] both checksums passed \n\n" ) ;
  }

  return corr ;
}
