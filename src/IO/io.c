/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "GLU_bswap.h"

// the question is ... Who checks the checksum?
void 
check_checksum( FILE *fprop )
{
  const int spinsize = NS * NS * NC * NC ;

  int site;
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;

  double *prop_buf = malloc( 2 * spinsize * sizeof(double) ) ;	

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
	   struct spinor *S )
{
  // size of a propagator at a site
  const int spinsize = NC * NC * NS * NS ;

  // read in site-by-site
  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( S[i].D , sizeof(double complex) , spinsize , fprop ) != 
	spinsize ) {
      printf( "Fread propagator failure \n" ) ;
      return FAILURE ;
    }
    // checksum that noise?
  }

  return SUCCESS ;
}

// Read NRQCD propagator time slice 
int 
read_nrprop( FILE *fprop, 
	     struct spinor *S )
{
  const int NR_NS = NS/2 ;
  const int spinsize = NC * NC * NR_NS * NR_NS ;

  //
  struct spinor P ;

  // read in site-by-site
  double *tmp = malloc( spinsize * 2 * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( tmp , sizeof(double complex), spinsize , fprop) != 
	spinsize ) {
      printf( "Fread propagator failure \n" ) ;
      free( tmp ) ;
      return FAILURE ;
    }

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
    // checksum that noise


    // change from nrel into chiral basis
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
    //  Note: There is another gamma_0 at the end!
    //

    for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {

	P.D[0][0].C[c1][c2] = 0.5 * (  S[i].D[0][0].C[c1][c2] + S[i].D[2][0].C[c1][c2] + S[i].D[0][2].C[c1][c2] + S[i].D[2][2].C[c1][c2] );
	P.D[1][0].C[c1][c2] = 0.5 * (  S[i].D[1][0].C[c1][c2] + S[i].D[3][0].C[c1][c2] + S[i].D[1][2].C[c1][c2] + S[i].D[3][2].C[c1][c2] );
	P.D[2][0].C[c1][c2] = 0.5 * ( -S[i].D[0][0].C[c1][c2] + S[i].D[2][0].C[c1][c2] - S[i].D[0][2].C[c1][c2] + S[i].D[2][2].C[c1][c2] );
	P.D[3][0].C[c1][c2] = 0.5 * ( -S[i].D[1][0].C[c1][c2] + S[i].D[3][0].C[c1][c2] - S[i].D[1][2].C[c1][c2] + S[i].D[3][2].C[c1][c2] );

	P.D[0][1].C[c1][c2] = 0.5 * (  S[i].D[0][1].C[c1][c2] + S[i].D[2][1].C[c1][c2] + S[i].D[0][3].C[c1][c2] + S[i].D[2][3].C[c1][c2] );
	P.D[1][1].C[c1][c2] = 0.5 * (  S[i].D[1][1].C[c1][c2] + S[i].D[3][1].C[c1][c2] + S[i].D[1][3].C[c1][c2] + S[i].D[3][3].C[c1][c2] );
	P.D[2][1].C[c1][c2] = 0.5 * ( -S[i].D[0][1].C[c1][c2] + S[i].D[2][1].C[c1][c2] - S[i].D[0][3].C[c1][c2] + S[i].D[2][3].C[c1][c2] );
	P.D[3][1].C[c1][c2] = 0.5 * ( -S[i].D[1][1].C[c1][c2] + S[i].D[3][1].C[c1][c2] - S[i].D[1][3].C[c1][c2] + S[i].D[3][3].C[c1][c2] );

	P.D[0][2].C[c1][c2] = 0.5 * ( -S[i].D[0][0].C[c1][c2] - S[i].D[2][0].C[c1][c2] + S[i].D[0][2].C[c1][c2] + S[i].D[2][2].C[c1][c2] );
	P.D[1][2].C[c1][c2] = 0.5 * ( -S[i].D[1][0].C[c1][c2] - S[i].D[3][0].C[c1][c2] + S[i].D[1][2].C[c1][c2] + S[i].D[3][2].C[c1][c2] );
	P.D[2][2].C[c1][c2] = 0.5 * (  S[i].D[0][0].C[c1][c2] - S[i].D[2][0].C[c1][c2] - S[i].D[0][2].C[c1][c2] + S[i].D[2][2].C[c1][c2] );
	P.D[3][2].C[c1][c2] = 0.5 * (  S[i].D[1][0].C[c1][c2] - S[i].D[3][0].C[c1][c2] - S[i].D[1][2].C[c1][c2] + S[i].D[3][2].C[c1][c2] );

	P.D[0][3].C[c1][c2] = 0.5 * ( -S[i].D[0][1].C[c1][c2] - S[i].D[2][1].C[c1][c2] + S[i].D[0][3].C[c1][c2] + S[i].D[2][3].C[c1][c2] );
	P.D[1][3].C[c1][c2] = 0.5 * ( -S[i].D[1][1].C[c1][c2] - S[i].D[3][1].C[c1][c2] + S[i].D[1][3].C[c1][c2] + S[i].D[3][3].C[c1][c2] );
	P.D[2][3].C[c1][c2] = 0.5 * (  S[i].D[0][1].C[c1][c2] - S[i].D[2][1].C[c1][c2] - S[i].D[0][3].C[c1][c2] + S[i].D[2][3].C[c1][c2] );
	P.D[3][3].C[c1][c2] = 0.5 * (  S[i].D[1][1].C[c1][c2] - S[i].D[3][1].C[c1][c2] - S[i].D[1][3].C[c1][c2] + S[i].D[3][3].C[c1][c2] );

      }
    }

    // put the result back into the spinor, also perform the final gamma_0
    /*for( c1 = 0 ; c1 < NC ; c1++ ) {
      for( c2 = 0 ; c2 < NC ; c2++ ) {

      S[i].D[0][0].C[c1][c2] = P.D[3][3].C[c1][c2];
      S[i].D[0][1].C[c1][c2] = P.D[3][2].C[c1][c2];
      S[i].D[0][2].C[c1][c2] = P.D[3][1].C[c1][c2];
      S[i].D[0][3].C[c1][c2] = P.D[3][0].C[c1][c2];

      S[i].D[1][0].C[c1][c2] = P.D[2][3].C[c1][c2];
      S[i].D[1][1].C[c1][c2] = P.D[2][2].C[c1][c2];
      S[i].D[1][2].C[c1][c2] = P.D[2][1].C[c1][c2];
      S[i].D[1][3].C[c1][c2] = P.D[2][0].C[c1][c2];

      S[i].D[2][0].C[c1][c2] = P.D[1][3].C[c1][c2];
      S[i].D[2][1].C[c1][c2] = P.D[1][2].C[c1][c2];
      S[i].D[2][2].C[c1][c2] = P.D[1][1].C[c1][c2];
      S[i].D[2][3].C[c1][c2] = P.D[1][0].C[c1][c2];

      S[i].D[3][0].C[c1][c2] = P.D[0][3].C[c1][c2];
      S[i].D[3][1].C[c1][c2] = P.D[0][2].C[c1][c2];
      S[i].D[3][2].C[c1][c2] = P.D[0][1].C[c1][c2];
      S[i].D[3][3].C[c1][c2] = P.D[0][0].C[c1][c2];

      }
      }*/

    // Re-shuffle source and sink indices
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    S[i].D[d1][d2].C[c1][c2] = P.D[d2][d1].C[c2][c1];
	  }
	}
      }
    }
  }
  free( tmp ) ;

  return SUCCESS ;
}




