/**
   @file mesonfile.c
   @brief little correlation file reader
 */
#include "common.h"

#include "correlators.h" // free_momcorr()
#include "reader.h"      // process_file()

// lattice geometry information
struct latt_info Latt ;

// 
enum { INFILE = 1 , OUTFILE = 2 } ;

// little code for accessing elements of our correlator files
int 
main( const int argc ,
      const char *argv[] )
{
  if( argc < 2 ) {
    return fprintf( stdout , "[TETAVG] Usage ./TETAVG {correlator file}"
		    "{outfile} \n" ) ;
  }

  // read the correlation file
  FILE *infile = fopen( argv[ INFILE ] , "rb" ) ;
  if( infile == NULL ) {
    fprintf( stderr , "[TETAVG] File %s does not exist\n" , argv[INFILE] ) ;
    return FAILURE ;
  }

  // set geometry to 0
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 0 ;
  }

  uint32_t NGSRC[1] = { 0 } , NGSNK[1] = { 0 } , NMOM[1] = { 0 } ;

  // allocate the momentum correlator
  struct mcorr **corr = NULL , **corravg ;
  struct veclist *momentum = NULL ;

  // read in the file and pack our structs
  corr = process_file( &momentum , infile , NGSRC , NGSNK , NMOM ) ;
  
  if( corr == NULL ) goto memfree ;

  corravg = allocate_momcorrs( (size_t)1 , NGSRC[0] , NMOM[0] ) ;

  // xyz average
  size_t GSRC , GSNK , p , t ;
  for( GSRC = 0 ; GSRC < NGSRC[0] ; GSRC++ ) {
    for( p = 0 ; p < NMOM[0] ; p++ ) {
      for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
	register double complex sum = 0.0 ;
	for( GSNK = 0 ; GSNK < NGSNK[0] ; GSNK++ ) {
	  sum += corr[GSRC][GSNK].mom[p].C[t] ; 
	}
	corravg[0][GSRC].mom[p].C[t] = sum / (double)NGSNK[0] ;
      }
    }
  }
  
  // write out the averaged list
  int nmom[1] = { (int)NMOM[0] } ;
  write_momcorr( argv[ OUTFILE ] , (const struct mcorr**)corravg ,
		 momentum , 1 , NGSRC[0] , nmom , "" ) ;
  
  free_momcorrs( corravg , 1 , NGSRC[0] , NMOM[0] ) ;
  
 memfree :

  // free the memory
  free_momcorrs( corr , NGSRC[0] , NGSNK[0] , NMOM[0] ) ;

  // free the momentum list
  free( momentum ) ;

  fclose( infile ) ;

  return SUCCESS ;
}
