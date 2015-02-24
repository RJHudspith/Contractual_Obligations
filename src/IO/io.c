/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "GLU_bswap.h"

// the question is ... Who checks the checksum?
<<<<<<< HEAD
int
check_checksum( FILE *fprop )
=======
void 
check_checksum( FILE *fprop , 
		const long int header )
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1
{
  // spinsize
  const int spinsize = NS * NS * NC * NC ;

  double *prop_buf = malloc( spinsize*2 * sizeof( double ) ) ;

<<<<<<< HEAD
  // crc accumulators
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;
=======
  double *prop_buf = malloc( 2 * spinsize * sizeof(double) ) ;	
  fseek(fprop, header, SEEK_SET) ;
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1

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
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( tmp , sizeof(double), spinsize*2 , fprop) != 
	spinsize*2 ) {
      printf( "Fread propagator failure \n" ) ;
      free( tmp ) ;
      return FAILURE ;
    }
<<<<<<< HEAD
=======

    // poke it into our struct
    int d1 , d2 , k = 0 ;
    for( d1 = 0 ; d1 < NS ; d1++ ) {
      for( d2 = 0 ; d2 < NS ; d2++ ) {
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
    // checksum that noise

>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1
  }
  free( tmp ) ;

  return SUCCESS ;
}



// Read NRQCD propagator time slice 
int 
read_nrprop( FILE *fprop, 
	   struct spinor *S ,
	   const long int header, 
	   const int tslice )
{

  int NR_NS=NS*0.5;
  const int spinsize = NC * NC * NR_NS * NR_NS ;
  long int jumper;
  struct spinor P;

  // Set position in file 
  jumper = header + VOL3 * tslice * spinsize *2 * sizeof(double);
  fseek(fprop, jumper, SEEK_SET);

  // read in site-by-site
  double *tmp = malloc( spinsize * 2 * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // Read in three timeslices from tslice 
    if( fread( tmp , sizeof(double), spinsize*2 , fprop) != 
	spinsize*2 ) {
      printf( "Fread propagator failure \n" ) ;
      free( tmp ) ;
      return FAILURE ;
    }

<<<<<<< HEAD
    int d1 , d2 , c1 , c2 , k = 0 ;
    // Now fill the non-zero components from the buffer
    for( d1 = 0 ; d1 < NR_NS ; d1++ ) {
      for( d2 = 0 ; d2 < NR_NS ; d2++ ) {
=======
	// Fill the buffer into the spinor
    int d1 , d2 , k = 0 ;
	int c1 , c2 ;
    // The NRQCD spinor has only 2x2 components
	// First, fill the spinor with zeros
    for( d1 = 0 ; d1 < NS ; d1++ ) { 
      for( d2 = 0 ; d2 < NS ; d2++ ) { 
		for( c1 = 0 ; c1 < NC ; c1++ ) {
	  	  for( c2 = 0 ; c2 < NC ; c2++ ) {
				S[i].D[d1][d2].C[c1][c2] = 0 + 0*I;
		  }
		} 	
	  }
	}

	// Now fill the non-zero components from the buffer
    for( d1 = 2 ; d1 < 2+NR_NS ; d1++ ) {
      for( d2 = 2 ; d2 < 2+NR_NS ; d2++ ) {
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1
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
<<<<<<< HEAD

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
=======
    // checksum that noise

  }
  free( tmp ) ;

  return SUCCESS ;
}
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1



// Read light propagator time slice and change to nrel basis
int 
read_chiral2nrel( FILE *fprop, 
	   struct spinor *S ,
	   const long int header, 
	   const int tslice )
{
  const int spinsize = NC * NC * NS * NS ;
  long int jumper;
  struct spinor P;
  int c1, c2;


  // Set position in file 
  jumper = header + VOL3 * tslice * spinsize *2 * sizeof(double);
  fseek(fprop, jumper, SEEK_SET);

  // read in site-by-site
  double *tmp = malloc( spinsize * 2 * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
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
	//int c1 , c2 ;
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    S[i].D[d1][d2].C[c1][c2] = tmp[ k ] + I * tmp[ k + 1 ] ; k+= 2 ;
	  }
	}
	#endif
      }
    }
    // checksum that noise

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
	// Re-shuffle source and sink indices
	for( d1 = 0 ; d1 < NS ; d1++ ) {
       for( d2 = 0 ; d2 < NS ; d2++ ) {
		  for( c1 = 0 ; c1 < NC ; c1++ ) {
       		for( c2 = 0 ; c2 < NC ; c2++ ) {
				
				S[i].D[d1][d2].C[c1][c2] = P.D[d1][d2].C[c1][c2];

			}
		  }
		}
	}
	

  }
  free( tmp ) ;

  return SUCCESS ;
}


