/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "GLU_bswap.h"

// the question is ... Who checks the checksum?
int
check_checksum( FILE *fprop )
{
  // spinsize
  const int spinsize = NS * NS * NC * NC ;

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
			(unsigned char*)prop_buf , 
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

// Read propagator time slice 
// should we accumulate the checksum? I think we should in fact I will probably do that next - J
int 
read_prop( FILE *fprop, 
	   struct spinor *S )
{
  const int spinsize = NC * NC * NS * NS ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( S[i].D , sizeof(double complex), spinsize , fprop) != 
	spinsize ) {
      printf( "Fread propagator failure \n" ) ;
      return FAILURE ;
    }
  }

  return SUCCESS ;
}

// Read NRQCD propagator time slice 
int 
read_nrprop( FILE *fprop, 
	     struct spinor *S )
{
  int NR_NS = NS/2 ;
  const int spinsize = NC * NC * NR_NS * NR_NS ;
  long int jumper;
  struct spinor P;

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

    int d1 , d2 , c1 , c2 , k = 0 ;
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

// Read light propagator time slice and change to nrel basis
int 
read_chiral2nrel( FILE *fprop, 
		  struct spinor *S )
{
  const int spinsize = NC * NC * NS * NS ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {

    struct spinor P ;
    // Read in three timeslices from tslice 
    if( fread( P.D , sizeof(double complex), spinsize , fprop ) != 
	spinsize ) {
      printf( "[C2NR] Fread propagator failure \n" ) ;
      return FAILURE ;
    }

    int c1, c2 , d1 , d2 ;

    // change from chiral into nrel basis
    //  
    //  chiral S = T * S * T^dagger
    // 
    //  here T:
    // 
    //        1/s    0       1/s     0
    //        0      1/s     0       1/s
    //        -1/s   0       1/s     0
    //        0      -1/s    0       1/s
    // 
    //  with 1/s = 1/sqrt(2)
    //       
    for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	S[i].D[0][0].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
	S[i].D[1][0].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
	S[i].D[2][0].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
	S[i].D[3][0].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

	S[i].D[0][1].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
	S[i].D[1][1].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
	S[i].D[2][1].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
	S[i].D[3][1].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );

	S[i].D[0][2].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
	S[i].D[1][2].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
	S[i].D[2][2].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
	S[i].D[3][2].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

	S[i].D[0][3].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
	S[i].D[1][3].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
	S[i].D[2][3].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
	S[i].D[3][3].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );	
      }
    }
    //
  }

  return SUCCESS ;
}


