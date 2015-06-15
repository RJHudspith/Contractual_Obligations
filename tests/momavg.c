/**
   @file mesonfile.c
   @brief little correlation file reader
 */
#include "common.h"

#include "corr_sort.h"    // generic sorting functions
#include "correlators.h"  // allocate and free corrs
#include "crc32.h"        // do the crc of the binary data
#include "cut_routines.h" // zero_veclist()
#include "geometry.h"     // init_geom()
#include "GLU_bswap.h"    // byte swaps if necessary
#include "GLU_timer.h"    // byte swaps if necessary

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

  // read the first on its own
  if( read_corr( corr[ 0 ][ 0 ].mom[ p ].C , cksuma , cksumb , 
		 infile , 0 ) == FAILURE ) {
    printf( "[IO] corr Read failure \n" ) ;
    return FAILURE ;
  }

  // read the next, checking we get the right LT
  size_t GSGK ;
  for( GSGK = 1 ; GSGK < NGSRC[0]*NGSNK[0] ; GSGK++ ) {
    const size_t GSRC = GSGK / NGSNK[0] ;
    const size_t GSNK = GSGK % NGSNK[0] ;
    // read the timeslice stuff
    uint32_t LT[1] ;
    if( FREAD32( LT , 1 , infile ) == FAILURE ) return FAILURE ;
    if( (int)LT[0] != L0 ) { 
      printf( "[IO] LT Read failure %d %d \n" , (int)LT[0] , L0 ) ; 
      return FAILURE ; 
    }
    // read the correlator
    if( read_corr( corr[ GSRC ][ GSNK ].mom[ p ].C , cksuma , cksumb , infile ,
		   GSNK + NGSNK[0] * GSRC ) == FAILURE ) {
      printf( "[IO] corr Read failure \n" ) ;
      return FAILURE ;
    }
  }

  return SUCCESS ;
}

// n^2 equivalents counter, data must be sorted
static int
count_equivalents( const struct veclist *list ,
		   const int NMOM )
{
  int m ;
  int Nequiv = 0 ;
  for( m = 0 ; m < NMOM ; m++ ) {
    int k ;
    for( k = 1 ; k < NMOM ; k++ ) { 
      if( ( m + k ) > ( NMOM  - 1 ) ) break ;
      if( list[m].nsq != list[k+m].nsq ) {
	break ;
      }
    }
    Nequiv ++ ;
    m += ( k - 1 ) ;
  }
  return Nequiv ;
}

static inline int
lt_nsq( const void *a , 
	const void *b )
{
  return (( *(const struct veclist*)a ).nsq 
	  < 
	  ( *(const struct veclist*)b ).nsq ) ;
}

static void
SORT( struct veclist *list ,
      struct mcorr **corr ,
      const int NMOM ,
      const int NGSRC ,
      const int NGSNK )
{
  // create a map
  size_t *map = malloc( NMOM * sizeof( size_t ) ) ;
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NMOM ; i++ ) {
    map[ i ] = i ;
  }

  merge_sort( list , map , sizeof( list[0] ) , sizeof( map[0] ), NMOM , lt_nsq ) ;

  // apply the sort to the correlators
  size_t GSGK ;
#pragma omp parallel for private(GSGK)
  for( GSGK = 0 ; GSGK < NGSRC * NGSNK ; GSGK++ ) {
    const size_t GSRC = GSGK / NGSNK ;
    const size_t GSNK = GSGK % NGSNK ;
    struct correlator *m = malloc( NMOM * sizeof( struct correlator ) ) ;
    size_t p ;
    for( p = 0 ; p < NMOM ; p++ ) {
      m[ p ].C = malloc( L0 * sizeof( double complex ) ) ;
      memcpy( m[ p ].C , corr[ GSRC ][ GSNK ].mom[ map[ p ] ].C , sizeof( double complex ) * L0 ) ;
    }
    for( p = 0 ; p < NMOM ; p++ ) {
      memcpy( corr[ GSRC ][ GSNK ].mom[ p ].C , m[ p ].C , sizeof( double complex ) * L0 ) ;
      free( m[ p ].C ) ;
    }
    free( m ) ;
  }		   
  free( map ) ;
  return ;
}

// equate two momentum correlators
static void
equate_momcorrs( struct mcorr **res ,
		 const struct mcorr **corr ,
		 const int momidx1 ,
		 const int momidx2 ,
		 const int NGSRC ,
		 const int NGSNK )
{
  int GSGK ;
#pragma omp parallel for private(GSGK)
  for( GSGK = 0 ; GSGK < ( NGSRC * NGSNK ) ; GSGK++ ) {
    const int GSRC = GSGK / NGSNK ;
    const int GSNK = GSGK % NGSNK ;
    memcpy( res[ GSRC ][ GSNK ].mom[ momidx1 ].C , 
	    corr[ GSRC ][ GSNK ].mom[ momidx2 ].C ,
	    L0 * sizeof( double complex ) ) ;
  }
}

// atomically add momentum correlators
static void
add_momcorrs( struct mcorr **res ,
	      const struct mcorr **corr ,
	      const int momidx1 ,
	      const int momidx2 ,
	      const int NGSRC ,
	      const int NGSNK )
{
  int GSGK ;
#pragma omp parallel for private(GSGK)
  for( GSGK = 0 ; GSGK < ( NGSRC * NGSNK ) ; GSGK++ ) {
    const int GSRC = GSGK / NGSNK ;
    const int GSNK = GSGK % NGSNK ;
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      res[ GSRC ][ GSNK ].mom[ momidx1 ].C[ t ] += 
	corr[ GSRC ][ GSNK ].mom[ momidx2 ].C[ t ] ;
    }
  }
}

// divide momentum correlator by a constant
static void
divide_constant( struct mcorr **res ,
		 const double constant ,
		 const int momidx1 ,
		 const int NGSRC ,
		 const int NGSNK )
{
  const double invconstant = 1.0 / constant ;
  int GSGK ;
#pragma omp parallel for private(GSGK)
  for( GSGK = 0 ; GSGK < ( NGSRC * NGSNK ) ; GSGK++ ) {
    const int GSRC = GSGK / NGSNK ;
    const int GSNK = GSGK % NGSNK ;
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      res[ GSRC ][ GSNK ].mom[ momidx1 ].C[ t ] *= invconstant ;
    }
  }
  return ;
}

struct mcorr**
momentum_average( struct veclist **avlist ,
		  int *Nequiv ,
		  const struct veclist *list ,
		  const struct mcorr **corr ,
		  const int NMOM ,
		  const int NGSRC ,
		  const int NGSNK )
{
  *Nequiv = count_equivalents( list , NMOM ) ;

  printf( "\n[AVE] %d equivalents \n" , *Nequiv ) ;

  // stop zero bytes malloc
  if( *Nequiv < 1 ) {
    return NULL ;
  }

  *avlist = malloc( *Nequiv * sizeof( struct veclist ) ) ;
  struct mcorr **corravg = allocate_momcorrs( NGSRC , NGSNK , *Nequiv ) ;

  int idx = 0 , m ;
  for( m = 0 ; m < NMOM ; m++ ) {

    // set the average list
    (*avlist)[ idx ] = list[ m ] ;
    
    // set average to corr
    equate_momcorrs( corravg , corr , idx , m , NGSRC , NGSNK ) ;

    #ifdef verbose
    printf( "Average :: ( %d %d ) [ %d %d %d ] -> [ %d %d %d ] :: ( %1.12e ) \n" , 
	    m , m ,
	    list[m].MOM[0] , list[m].MOM[1] , list[m].MOM[2] , 
	    list[m].MOM[0] , list[m].MOM[1] , list[m].MOM[2] , 
	    creal( corr[5][5].mom[m].C[1] ) ) ;
    #endif

    // loop through equivalent momenta again
    int k , nequiv = 1 ;
    for( k = 1 ; k < NMOM ; k++ ) { 
      if( ( m + k ) > ( NMOM  - 1 ) ) break ;
      if( list[m].nsq != list[k+m].nsq ) {
	break ;
      } else {
	#ifdef verbose
	printf( "Average :: ( %d, %d ) [ %d %d %d ] -> [ %d %d %d ] :: ( %1.12e ) \n" , 
		m , k+m ,
		list[m].MOM[0] , list[m].MOM[1] , list[m].MOM[2] , 
		list[m+k].MOM[0] , list[m+k].MOM[1] , list[m+k].MOM[2] , 
		creal( corr[5][5].mom[k+m].C[1] ) ) ;
	#endif
	add_momcorrs( corravg , corr , idx , k+m , NGSRC , NGSNK ) ;
	nequiv++ ;
      }
    }
    divide_constant( corravg , (double)(nequiv) , idx , NGSRC , NGSNK ) ;

    #ifdef verbose
    printf( "AVERAGE :: %d %1.12e \n" , idx , creal( corravg[5][5].mom[idx].C[1] ) ) ;
    #endif

    idx ++ ;
    m += ( k - 1 ) ;
  }

  return corravg ;
}

// print to stdout the momentum correlators
#ifdef verbose
static void
print_mcorrs( const struct veclist *momentum , 
	      const struct mcorr **corr ,
	      const int NMOM ,
	      const int GSRC ,
	      const int GSNK )
{
  int p ;
  for( p = 0 ; p < NMOM ; p++ ) {
    int mu ;
    printf( "(" ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      printf( " %d " , momentum[p].MOM[mu] ) ;
    }
    printf( ") %d \n" , momentum[p].nsq ) ;
    int t ;
    for( t = 0 ; t < L0 ; t++ ) {
      printf( "%d %e %e \n" , t , 
	      creal( corr[ GSRC ][ GSNK ].mom[p].C[t] ) ,
	      cimag( corr[ GSRC ][ GSNK ].mom[p].C[t] ) ) ;
    }
    printf( "\n" ) ;
  }
  return ;
}
#endif

// write out a whole heap of averaged correlators
static void
write_averages( const struct veclist *list ,
		const struct mcorr **corr ,
		const char *outname ,
		const int NMOM , 
		const int NGSRC ,
		const int NGSNK )
{
  // create a dummy correlator
  int *dNMOM = malloc( sizeof( int ) ) ;
  struct veclist *dlist = (struct veclist*)zero_veclist( dNMOM , ND-1 , GLU_FALSE ) ;
  struct mcorr **dcorr = allocate_momcorrs( NGSRC , NGSNK , dNMOM[0] ) ;

  // loop averaged momenta
  int p ;
  for( p = 0 ; p < NMOM ; p++ ) {

    // set the dummy correlator such that dcorr[0] = corr[p]
    equate_momcorrs( dcorr , corr , 0 , p , NGSRC , NGSNK ) ;

    // set the dlist
    dlist[0] = list[p] ;

    // write out sequential "zero-mom" files 
    char str[ 256 ] ;
    sprintf( str , "%s.%d.bin" , outname , p ) ;
    write_momcorr( str , (const struct mcorr **)dcorr , dlist , NGSRC , NGSNK , dNMOM ) ;
  }

  // free the dummies
  free_momcorrs( dcorr , NGSRC , NGSNK , dNMOM[0] ) ;
  free( dlist ) ;
  free( dNMOM ) ;

  return ;
}

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc != 3 ) {
    return printf( "usage ./MESONS {correlator file} {outfile}\n" ) ;
  }
  printf( "[IO] Outputting to %s \n" , argv[2] ) ;

  // read the correlation file
  FILE *infile = fopen( argv[1] , "rb" ) ;
  if( infile == NULL ) {
    printf( "File %s does not exist\n" , argv[1] ) ;
    return -1 ;
  }

  start_timer( ) ;

  uint32_t magic[1] , NGSRC[1] = { 0 } ;
  uint32_t NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

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

  struct veclist *momentum = malloc( NMOM[0] * sizeof( struct veclist ) ) ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL , **corravg = NULL ;

  // momentum averages
  int Nequiv = 0 ;
  struct veclist *avlist = NULL ;

  // checksums
  uint32_t cksuma = 0 , cksumb = 0 , csum[ 2 ] = { 0 , 0 } ; 

  int p ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    uint32_t n[ ND ] ;
    if( FREAD32( n , ND , infile ) == FAILURE ) goto memfree ;
    if( n[ 0 ] != ND-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , ND-1 ) ;
      goto memfree ;
    }
    int mu ;
    momentum[ p ].nsq = 0 ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      momentum[ p ].MOM[ mu ] = (int)n[ 1 + mu ] ;
      momentum[ p ].nsq += ( n[ 1 + mu ] * n[ 1 + mu ] ) ;
    }
    momentum[ p ].idx = p ;
  }

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

  print_time( ) ;

  // check our checksums
  if( FREAD32( csum , 2 , infile ) == FAILURE ) goto memfree ;
  if( csum[0] != cksuma || csum[1] != cksumb ) {
    printf( "Mismatched checksums ! %x %x %x %x\n" , csum[0] , csum[1] , cksuma , cksumb ) ;
    goto memfree ;
  } 

  printf( "\n[CHECKSUM] both checksums passed \n" ) ;

  print_time( ) ;

  // sort the momenta && mom correlator
  SORT( momentum , corr , NMOM[0] , NGSRC[0] , NGSNK[0] ) ;

  printf( "\n[SORT] momentum list sorted \n" ) ;

  print_time( ) ;

  // average
  corravg = momentum_average( &avlist , &Nequiv ,
			      momentum , (const struct mcorr**)corr ,
			      NMOM[0] , NGSRC[0] , NGSNK[0] ) ;

  printf( "\n[AVE] correlators averaged \n" ) ;

  print_time( ) ;

  if( corravg == NULL ) goto memfree ;

  // have a look at some of the correlators
#ifdef verbose
  print_mcorrs( avlist , (const struct mcorr**)corravg , Nequiv , 5 , 5 ) ;
#endif

  // split the averaged results into separate files for ease of reading
  write_averages( avlist , (const struct mcorr**)corravg , argv[2] , 
		  Nequiv , NGSRC[0] , NGSNK[0] ) ;

  printf( "\n[MOMAVG] all finished \n" ) ;

  print_time( ) ;

 memfree :

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  // free the average momentum list
  free( avlist ) ;

  // free the average correlator
  free_momcorrs( corravg , NGSRC[0] , NGSNK[0] , Nequiv ) ;

  // close the file
  fclose( infile ) ;

  return SUCCESS ;
}

