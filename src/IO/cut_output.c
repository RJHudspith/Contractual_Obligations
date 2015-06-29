/**
   @file cut_output.c
   @brief outputs our cut data
 */
#include "common.h"

#include "crc32c.h"        // crc 

// Write the gluonic two point function to the file Ap
static void
write_binary( FILE *__restrict Ap , 
	      const double *__restrict g2 , 
	      const int NMOM[ 1 ] ) 
{
  fwrite( NMOM , sizeof(uint32_t) , 1 , Ap ) ;
  //fwrite( g2 , sizeof(double) , NMOM[0] , Ap ) ;
  // checksum it
  uint32_t cksuma = 0 , cksumb = 0 ;
  int i ;
  for( i = 0 ; i < NMOM[0] ; i++ ) {
    fwrite( g2 + i , sizeof(double) , 1 , Ap ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , i , 
			       (void*)( g2 + i ) , sizeof( double  ) ) ;
  }
  const uint32_t csum[ 2 ] = { cksuma , cksumb } ;
  fwrite( csum , sizeof( uint32_t ) , 2 , Ap ) ;

  return ;
}

// Support for writing our momentum list
void
write_mom_veclist( FILE *__restrict Ap , 
		   const int *__restrict num_mom , 
		   const struct veclist *__restrict list ,
		   const int DIR )
{
  const int stride = ( DIR + 1 ) * num_mom[ 0 ] ;
  int *kall = (int*)malloc( stride * sizeof(int) ) ;

  fwrite( num_mom , sizeof(int) , 1 , Ap ) ;

  int i ;
  //write momenta
#pragma omp parallel for private(i)
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    int n = i * ( DIR + 1 ) ;
    kall[ n ] = DIR ; 
    n++ ; 
    int a ;
    for( a = 0  ;  a < DIR  ;  a++ ) {
      kall[ n ] = list[i].MOM[a] ; 
      n++ ; 
    }
  }
  fwrite( kall , sizeof(int) , stride , Ap ) ;
  free( kall ) ;
  return ;
}

//
void
write_momspace_data( const char *filename ,  
		     const int *__restrict NMOM ,
		     const double *__restrict data ,
		     const struct veclist *__restrict list ,
		     const int DIR )
{
  printf( "[IO] writing momspace data to %s \n" , filename ) ;

  FILE *outfile = fopen( filename , "wb" ) ;

  // write a magic number
  uint32_t magic[ 1 ] = { VPF_MAGIC } ;
  fwrite( magic , sizeof( uint32_t ) , 1 , outfile ) ;

  // writes the momentum list
  // in terms of Fourier modes px py pz pt
  write_mom_veclist( outfile , NMOM , list , DIR ) ;

  // writes out the binary data
  write_binary( outfile , data , NMOM ) ;

  fclose( outfile ) ;

  return ;
}

