/**
   @file momavg.c
   @brief reads correlation file and averages equivalent fourier modes
 */
#include "common.h"

#include "corr_sort.h"    // generic sorting functions
#include "correlators.h"  // free_momcorrs()
#include "cut_routines.h" // zero_veclist()
#include "GLU_timer.h"    // print_time()
#include "reader.h"       // process_file()

// lattice geometry
struct latt_info Latt ;

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

// less than veclist sort
static inline int
lt_nsq( const void *a , 
	const void *b )
{
  return ( ( *(const struct veclist*)a ).nsq < 
	   ( *(const struct veclist*)b ).nsq ) ;
}

// sort the data by the veclist
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
      m[ p ].C = malloc( LT * sizeof( double complex ) ) ;
      memcpy( m[ p ].C , corr[ GSRC ][ GSNK ].mom[ map[ p ] ].C , sizeof( double complex ) * LT ) ;
    }
    for( p = 0 ; p < NMOM ; p++ ) {
      memcpy( corr[ GSRC ][ GSNK ].mom[ p ].C , m[ p ].C , sizeof( double complex ) * LT ) ;
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
	    LT * sizeof( double complex ) ) ;
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
    for( t = 0 ; t < LT ; t++ ) {
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
    for( t = 0 ; t < LT ; t++ ) {
      res[ GSRC ][ GSNK ].mom[ momidx1 ].C[ t ] *= invconstant ;
    }
  }
  return ;
}

// perform n^2 momentum average
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

  fprintf( stdout , "\n[AVE] %d equivalents \n" , *Nequiv ) ;

  // stop zero bytes malloc
  if( *Nequiv < 1 ) return NULL ;

  *avlist = malloc( *Nequiv * sizeof( struct veclist ) ) ;
  struct mcorr **corravg = allocate_momcorrs( NGSRC , NGSNK , *Nequiv ) ;

  int idx = 0 , m ;
  for( m = 0 ; m < NMOM ; m++ ) {

    // set the average list
    (*avlist)[ idx ] = list[ m ] ;
    
    // set average to corr
    equate_momcorrs( corravg , corr , idx , m , NGSRC , NGSNK ) ;

    #ifdef verbose
    fprintf( stdout , "Average :: ( %d %d ) [ %d %d %d ] "
	     "-> [ %d %d %d ] :: ( %1.12e ) \n" , 
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
	fprintf( stdout , "Average :: ( %d, %d ) [ %d %d %d ] "
		 "-> [ %d %d %d ] :: ( %1.12e ) \n" , 
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
    fprintf( stdout , "AVERAGE :: %d %1.12e \n" , idx , 
	     creal( corravg[5][5].mom[idx].C[1] ) ) ;
    #endif

    idx ++ ;
    m += ( k - 1 ) ;
  }

  return corravg ;
}

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
    write_momcorr( str , (const struct mcorr **)dcorr , dlist , 
		   NGSRC , NGSNK , dNMOM , "" ) ;
  }

  // free the dummies
  free_momcorrs( dcorr , NGSRC , NGSNK , dNMOM[0] ) ;
  free( dlist ) ;
  free( dNMOM ) ;

  return ;
}

// enum for the IO
enum { INFILE = 1 , SPLIT = 2 , OUTFILE = 3 } ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  // check the number of command line arguments
  if( argc != 4 ) {
    return fprintf( stdout , "[MOMAVG] "
		    "usage ./MESONS {correlator file} {split} {outfile}\n" 
		    "{split} :: true/false\n" ) ;
  }
  fprintf( stdout , "[IO] Outputting to %s \n" , argv[ OUTFILE ] ) ;

  // check in and out are not the same
  if( are_equal( argv[ INFILE ] , argv[ OUTFILE ] ) ) {
    return FAILURE ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[ INFILE ] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[IO] File %s does not exist\n" , argv[ INFILE ] ) ;
    return FAILURE ;
  }

  // are we splitting up the files?
  GLU_bool split = GLU_FALSE ;
  if( are_equal( argv[ SPLIT ] , "true" ) ) {
    split = GLU_TRUE ;
  }

  start_timer( ) ;

  struct veclist *momentum = NULL ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL , **corravg = NULL ;

  // momentum averages
  int Nequiv = 0 ;
  struct veclist *avlist = NULL ;

  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;

  fprintf( stdout , "[IO] (N_SRC, N_SNK, NMOM) :: (%d, %d, %d)\n" , 
	   NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  fprintf( stdout , "\n[IO] file read \n" ) ;

  print_time( ) ;

  // sort the momenta && mom correlator
  if( NMOM[0] > 0 ) {
    SORT( momentum , corr , NMOM[0] , NGSRC[0] , NGSNK[0] ) ;
  } else {
    goto memfree ;
  }

  fprintf( stdout , "\n[SORT] momentum list sorted \n" ) ;

  print_time( ) ;

  // average
  corravg = momentum_average( &avlist , &Nequiv ,
			      momentum , (const struct mcorr**)corr ,
			      NMOM[0] , NGSRC[0] , NGSNK[0] ) ;

  fprintf( stdout , "\n[AVE] correlators averaged \n" ) ;

  print_time( ) ;

  if( corravg == NULL ) goto memfree ;

  if( split == GLU_TRUE ) {
    // split the averaged results into separate files for ease of reading
    write_averages( avlist , (const struct mcorr**)corravg , argv[3] , 
		    Nequiv , NGSRC[0] , NGSNK[0] ) ;
  } else {
    // write into a file
    int NMOM[1] = { Nequiv } ;
    write_momcorr( argv[ OUTFILE ] , (const struct mcorr **)corravg , 
		   avlist , NGSRC[0] , NGSNK[0] , (const int*)NMOM , "" ) ;
  }

  fprintf( stdout , "\n[MOMAVG] all finished \n" ) ;

  print_time( ) ;

 memfree :

  // free the memory
  if( NMOM[0] > 0 ) {
    free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;
  }

  // free the momentum list
  free( momentum ) ;

  // free the average momentum list
  free( avlist ) ;

  // free the average correlator
  if( Nequiv > 0 ) {
    free_momcorrs( corravg , NGSRC[0] , NGSNK[0] , Nequiv ) ;
  }

  // close the file
  fclose( infile ) ;

  return SUCCESS ;
}

