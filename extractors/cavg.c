/**
   @file cavg.c
   @brief average correlators
 */
#include "common.h"

#include "correlators.h"  // free_momcorrs()
#include "cut_routines.h" // zero_veclist()
#include "GLU_timer.h"    // print_time()
#include "reader.h"       // process_file()

// struct for storing our correlator files
struct correlator_file {
  struct mcorr **corr ;
  struct veclist *momentum ;
} ;

// lattice geometry
struct latt_info Latt ;

// little code for averaging correlation functions
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 3 ) {
    return fprintf( stdout , "[CAVG] "
		    "usage ./CAVG {c1} ... {c_n} {avg}\n" ) ;
  }
  fprintf( stdout , "[CAVG] Outputting to %s \n" , argv[ argc-1 ] ) ;

  const size_t Navg = (const size_t)argc-2 ;
  fprintf( stdout , "[CAVG] averaging %d files \n" , argc-2 ) ;

  size_t n ;
  struct correlator_file *cfiles = malloc( (Navg) * sizeof( struct correlator_file ) ) ;

  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  FILE *infile = fopen( argv[ 1 ] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[IO] File %s does not exist\n" , argv[ 1 ] ) ;
    goto memfree ;
  }
  cfiles[0].corr = process_file( &cfiles[0].momentum ,
				 infile , NGSRC , NGSNK , NMOM ) ;
  fprintf( stdout , "[IO] (N_SRC, N_SNK, NMOM) :: (%d, %d, %d)\n" , 
	   NGSRC[0] , NGSNK[0] , NMOM[0] ) ;
  fprintf( stdout , "\n[IO] file read \n" ) ;
  fclose( infile ) ;
  
  for( n = 1 ; n < Navg ; n++ ) {
    FILE *tinfile = fopen( argv[ n+1 ] , "rb" ) ;
    if( tinfile == NULL ) {
      fprintf( stderr , "[IO] File %s does not exist\n" , argv[ n+1 ] ) ;
      goto memfree ;
    }
    uint32_t tNGSRC[1] = { 0 } , tNGSNK[1] = { 0 } , tNMOM[1] = { 0 } ;
    cfiles[n].corr = process_file( &cfiles[n].momentum ,
				   tinfile , tNGSRC , tNGSNK , tNMOM ) ;
    fprintf( stdout , "[IO] (N_SRC, N_SNK, NMOM) :: (%d, %d, %d)\n" , 
	     tNGSRC[0] , tNGSNK[0] , tNMOM[0] ) ;
    fprintf( stdout , "\n[IO] file read \n" ) ;
    fclose( tinfile ) ;
    if( tNGSRC[0] != NGSRC[0] ||
	tNGSNK[0] != NGSNK[0] ||
	tNMOM[0] != NMOM[0] ) {
      fprintf( stderr , "[IO] files have different SRC/SNK/NMOM\n" ) ;
      goto memfree ;
    }
  }

  // do the average into element 0
  size_t GSGK ;
  #pragma omp parallel for private(GSGK)
  for( GSGK = 0 ; GSGK < ( NGSRC[0] * NGSNK[0] ) ; GSGK++ ) {
    const size_t GSRC = GSGK / NGSNK[0] ;
    const size_t GSNK = GSGK % NGSNK[0] ;
    size_t t , n , p ;
    for( p = 0 ; p < NMOM[0] ; p++ ) {
      for( t = 0 ; t < LT ; t++ ) {
	for( n = 1 ; n < Navg ; n++ ) {
	  cfiles[0].corr[ GSRC ][ GSNK ].mom[p].C[ t ] +=
	    cfiles[n].corr[ GSRC ][ GSNK ].mom[p].C[ t ] ;
	}
	cfiles[0].corr[ GSRC ][ GSNK ].mom[p].C[ t ] /= Navg ;
      }
    }
  }

  double twist_zero[ ND ] ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    twist_zero[ mu ] = 0.0 ;
  }
  // write the averaged correlator
  write_momcorr( argv[ argc-1 ] , (const struct mcorr**)cfiles[0].corr , 
		 cfiles[0].momentum , twist_zero , NGSRC[0] , NGSNK[0] ,
		 (const int*)NMOM , "" ) ;
  
 memfree :
  
  for( n = 0 ; n < Navg ; n++ ) {
    if( NMOM[0] > 0 ) {
      if( cfiles[n].corr != NULL ) {
	free_momcorrs( cfiles[n].corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;
      }
    }
    if( cfiles[n].momentum != NULL ) {
      free( cfiles[n].momentum ) ;
    }
  }
  
  return SUCCESS ;
}
