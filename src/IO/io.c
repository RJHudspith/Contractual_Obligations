/**
   @file io.c
   @file readers and alike

   TODO checksum accumulator in prop_read
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
read_prop( FILE *__restrict fprop, 
	   struct spinor *__restrict S ,
	   const long int header, 
	   const int tslice )
{
  const int spinsize = NC * NC * NS * NS ;
  long int jumper;

  // Set position in file 
  jumper = header + VOL3 * tslice * spinsize *2 * sizeof(double);
  fseek(fprop, jumper, SEEK_SET);

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // our struct is bytewise equivalent in C-order so we should be able 
    // to read it all in directly
    if( fread( S[i].D , sizeof(double complex), spinsize , fprop) != 
	spinsize ) {
      printf( "Fread propagator failure \n" ) ;
      return FAILURE ;
    }
    // accumulate the checksum 
  }

  return SUCCESS ;
}

