/**
   @file pent_avg.c
   @brief averages the forward solution with the opposite correctly projected parity backward one
 */
#include "common.h"

#include "correlators.h" // free_momcorr()
#include "reader.h"      // process_file()

// lattice geometry information
struct latt_info Latt ;

// 
enum { INFILE = 1 , TSRC = 2 , OUTFILE = 3 } ;

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc != 4 ) {
    return fprintf( stdout , "[PENTAVG] Usage ./PENTAVG {correlator file} tsrc"
		    " {outfile} \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[ INFILE ] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[PENTAVG] File %s does not exist\n" , argv[INFILE] ) ;
    return FAILURE ;
  }

  // set geometry to 0
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 0 ;
  }

  const int tsrc = atoi( argv[ TSRC ] ) ;
  
  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL , **corravg ;
  struct veclist *momentum = NULL ;

  // read in the file and pack our structs
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;
  
  if( corr == NULL ) goto memfree ;

  corravg = allocate_momcorrs( (size_t)1 , NGSNK[0] , NMOM[0] ) ;

  // xyz average
  size_t GSNK , p , t ;
  for( GSNK = 0 ; GSNK < NGSNK[0] ; GSNK++ ) {
    for( p = 0 ; p < NMOM[0] ; p++ ) {
      //
      corravg[0][GSNK].mom[p].C[0] =
	0.5*( corr[0][GSNK].mom[p].C[ 0 ] + corr[1][GSNK].mom[p].C[ 0 ] ) ;
      
      
      for( t = 1 ; t < LT ; t++ ) {
	double complex t1 , t2 ;
	
	if( t <= tsrc ) {
	  t1 =  corr[0][GSNK].mom[p].C[ LT-t ] ;
	} else {
	  t1 = -corr[0][GSNK].mom[p].C[ LT-t ] ;
	}
	//
	if( t >= (LT - tsrc) ) {
	  t2 = -corr[1][GSNK].mom[p].C[t] ;
	} else {
	  t2 =  corr[1][GSNK].mom[p].C[t] ;
	}
	corravg[0][GSNK].mom[p].C[t] = 0.5*(t1+t2) ;
      }
      //
    }
  }
  
  // write out the averaged list
  int nmom[1] = { (int)NMOM[0] } ;
  write_momcorr( argv[ OUTFILE ] , (const struct mcorr**)corravg ,
		 momentum , NULL , 1 , NGSNK[0] , nmom , "" ) ;
  free_momcorrs( corravg , 1 , NGSNK[0] , NMOM[0] ) ;
  
 memfree :

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  fclose( infile ) ;

  return SUCCESS ;
}
