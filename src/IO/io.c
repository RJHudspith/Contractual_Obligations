/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "crc32.h"
#include "GLU_bswap.h"

// the question is ... Who checks the checksum?
int
check_checksum( FILE *fprop )
{
  // spinsize
  const int spinsize = NSNS * NCNC ;

  double *prop_buf = malloc( spinsize*2 * sizeof( double ) ) ;

  // crc accumulators
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;

  // ok so this must be a geometry/dirac thing
  int site;
  for( site = 0 ; site < LVOLUME ; site++ ) {
    if( fread( prop_buf , sizeof(double) , spinsize*2 , fprop ) != spinsize*2 ) {
      printf( "[IO] fread failure ... exiting \n" ) ;
      return FAILURE ;
    }
    // we need to know what end is up, this is actually kinda tricky
    DML_checksum_accum( &CRCsum29 , &CRCsum31 , site , 
			(char*)prop_buf , 
			2 * spinsize * sizeof(double) ) ;
  }
  
  // we are at the end of the file now so we can fscanf for the value I will have to grok the output
  uint32_t rCRCsum29 , rCRCsum31 ;
  if( fscanf( fprop , "%x %x" , &rCRCsum29 , &rCRCsum31 ) != 2 ) {
    printf( "[IO] file read failure at the checksums \n" ) ;
    free( prop_buf ) ;
    return FAILURE ;
  }

  if( CRCsum29 != rCRCsum29 || CRCsum31 != rCRCsum31 ) {
    printf( "[IO] mismatched checksums \n" ) ;
    printf( "[IO] Computed Checksums %x %x\n" , CRCsum29 , CRCsum31 ) ;
    printf( "[IO] File Read Checksums %x %x\n" , rCRCsum29 , rCRCsum31 ) ;
    free( prop_buf ) ;
    return FAILURE ;
  }

  return SUCCESS ;
}

// Read light propagator on a time slice 
// should we accumulate the checksum? Probably
static int 
read_chiralprop( FILE *fprop, 
		 struct spinor *S )
{
  const int spinsize = NCNC * NSNS ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in propagator on a timeslice
    if( fread( S[i].D , sizeof(double complex), spinsize , fprop) != 
	spinsize ) {
      printf( "Fread propagator failure \n" ) ;
      return FAILURE ;
    }
  }

  return SUCCESS ;
}

// Read NRQCD propagator time slice 
static int 
read_nrprop( FILE *fprop, 
	     struct spinor *S )
{
  const int NR_NS = NS >> 1 ;
  const int spinsize = NCNC * NR_NS * NR_NS ;

  // read in site-by-site
  double *tmp = malloc( spinsize * 2 * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( tmp , sizeof(double complex), spinsize , fprop) != 
	spinsize ) {
      printf( "[NR] Fread propagator failure \n" ) ;
      free( tmp ) ;
      return FAILURE ;
    }

    int d1 , d2 , k = 0 ;
    // Now fill the non-zero components from the buffer
    for( d1 = 0 ; d1 < NR_NS ; d1++ ) {
      for( d2 = 0 ; d2 < NR_NS ; d2++ ) {
	#if NC == 3
	S[i].D[d1][d2].C[0][0] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[0][1] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[0][2] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[1][0] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[1][1] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[1][2] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[2][0] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[2][1] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	S[i].D[d1][d2].C[2][2] = tmp[k] + I * tmp[k+1] ; k+=2 ;
	#else
	int c1 , c2 ;
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    S[i].D[d1][d2].C[c1][c2] = tmp[ k ] + I * tmp[ k + 1 ] ; k+= 2 ;
	  }
	}
	#endif
      }
    }
  }
  free( tmp ) ;

  return SUCCESS ;
}

// thin wrapper for propagator reading
int
read_prop( struct propagator prop ,
	   struct spinor *S )
{
  switch( prop.basis ) {
  case CHIRAL :
    return read_chiralprop( prop.file , S ) ;
  case NREL :
    return read_nrprop( prop.file , S ) ;
  }
  return FAILURE ;
}

