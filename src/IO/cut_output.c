/**
   @file cut_output.c
   @brief outputs our cut data
 */
#include "common.h"

#include "crc32c.h"        // crc 

// Write the gluonic two point function to the file Ap
static void
write_binary( FILE *Ap , 
	      const double *g2 , 
	      const int NMOM[ 1 ] ) 
{
  fwrite( NMOM , sizeof(uint32_t) , 1 , Ap ) ;
  
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
write_mom_veclist( FILE *Ap , 
		   const int *num_mom , 
		   const struct veclist *list ,
		   const int DIR )
{
  size_t i ;
  fwrite( num_mom , sizeof(int) , 1 , Ap ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    fwrite( &DIR , sizeof(int) , 1 , Ap ) ;
    fwrite( list[i].MOM , sizeof(double) , DIR , Ap ) ;
  }
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
  fprintf( stdout , "[IO] writing momspace data to %s \n" , filename ) ;

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

