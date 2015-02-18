/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "GLU_bswap.h"

// the question is ... Who checks the checksum?
void 
check_checksum( FILE *fprop , 
		const long int header )
{
  const int spinsize = NS * NS * NC * NC ;

  int site;
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;

  double *prop_buf = malloc( 2 * spinsize * sizeof(double) ) ;	
  fseek(fprop, header, SEEK_SET) ;

  for( site = 0 ; site < VOL4 ; site++ ) {
    if( fread(prop_buf, sizeof(double), spinsize*2 , fprop) != spinsize*2 ) {
      printf( "fread failure ... exiting \n" ) ;
      exit( -1 ) ;
    }
    // we need to know what end is up, this is actually kinda tricky
    bswap_64( 2*spinsize , prop_buf ) ; 
    DML_checksum_accum( &CRCsum29, &CRCsum31, site, (unsigned char*)(prop_buf), 2*spinsize*sizeof(double) );
  }	
  
  // computed
  printf("[Computed Checksums] %x %x\n", CRCsum29 , CRCsum31 ) ;
  
  /*
  // we are at the end of the file now so we can fscanf for the value I will have to grok the output
  uint32_t rCRCsum29 , rCRCsum31 ;
  fscanf( file , "%x %x" , &rCRCsum29 , &rCRCsum31 ) ;
  printf("[File Read Checksums] %x %x\n", CRCsum29 , CRCsum31 ) ;
  */

  free( prop_buf ) ;

  return ;
}

// Read propagator time slice 
// should we accumulate the checksum? I think we should in fact I will probably do that next - J
int 
read_prop( FILE *fprop, 
	   struct spinor *S ,
	   const long int header, 
	   const int tslice )
{
  const int spinsize = NC * NC * NS * NS ;
  long int jumper;

  // Set position in file 
  jumper = header + VOL3 * tslice * spinsize *2 * sizeof(double);
  fseek(fprop, jumper, SEEK_SET);

  // read in site-by-site
  double *tmp = malloc( spinsize * 2 * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < VOL3 ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( tmp , sizeof(double), spinsize*2 , fprop) != 
	spinsize*2 ) {
      printf( "Fread propagator failure \n" ) ;
      free( tmp ) ;
      return FAILURE ;
    }
    // poke it into our struct
    int d1 , d2 , k = 0 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	#if NC == 3
	S[i].D[d1][d2].C[0][0] = tmp[k++] ;
	S[i].D[d1][d2].C[0][1] = tmp[k++] ;
	S[i].D[d1][d2].C[0][2] = tmp[k++] ;
	S[i].D[d1][d2].C[1][0] = tmp[k++] ;
	S[i].D[d1][d2].C[1][1] = tmp[k++] ;
	S[i].D[d1][d2].C[1][2] = tmp[k++] ;
	S[i].D[d1][d2].C[2][0] = tmp[k++] ;
	S[i].D[d1][d2].C[2][1] = tmp[k++] ;
	S[i].D[d1][d2].C[2][2] = tmp[k++] ;
	#else
	int c1 , c2 ;
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    S[i].D[d1][d2].C[c1][c2] = tmp[ k++ ] ;
	  }
	}
	#endif
      }
    }
    // checksum that noise

  }
  free( tmp ) ;

  return SUCCESS ;
}

