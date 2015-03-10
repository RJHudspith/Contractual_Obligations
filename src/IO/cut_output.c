/**
   @file cut_output.c
   @brief outputs our cut data
 */

#include "common.h"

// Write the gluonic two point function to the file Ap
static void
write_binary( FILE *__restrict Ap , 
	      const double *__restrict g2 , 
	      const int NMOM[ 1 ] ) 
{
  fwrite( NMOM , sizeof(uint32_t) , 1 , Ap ) ;
  fwrite( g2 , sizeof(double) , NMOM[0] , Ap ) ;
  return ;
}

// Support for writing our momentum list
static void
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

  // writes the momentum list
  // in terms of Fourier modes px py pz pt
  write_mom_veclist( outfile , NMOM , list , DIR ) ;

  // writes out the binary data
  write_binary( outfile , data , NMOM ) ;

  fclose( outfile ) ;
  return ;
}

// file writer
void
write_tmoments( const double **tcorr ,
		const char *filename ,
		const current_type current ,
		const vector_axial VA )
{
  char str[ 256 ] ;
  switch( current ) {
  case CONSERVED_LOCAL :
    switch( VA ) {
    case AXIAL :
      sprintf( str , "%s.CALA.tcorr.bin" , filename ) ;
      break ;
    case VECTOR :
      sprintf( str , "%s.CVLV.tcorr.bin" , filename ) ;
      break ;
    } break ;
  case LOCAL_LOCAL :
    switch( VA ) {
    case AXIAL :
      sprintf( str , "%s.LALA.tcorr.bin" , filename ) ;
      break ;
    case VECTOR :
      sprintf( str , "%s.LVLV.tcorr.bin" , filename ) ;
      break ;
    } break ;
  }
  printf( "[IO] writing temporal correlators to %s \n" , str ) ;

  FILE *outfile = fopen( str , "wb" ) ;

  // write a magic number
  const uint32_t magic[ 1 ] = { VPF_MAGIC } ;

  // write out the magic number
  fwrite( magic , sizeof( uint32_t ) , 1 , outfile ) ;

  // write out ND*ND
  uint32_t in[ 1 ] = { ND * ND } ;
  fwrite( in , sizeof( uint32_t ) , 1 , outfile ) ;

  // write out each munu loop
  int munu ;
  for( munu = 0 ; munu < ND*ND ; munu++ ) {
    in[ 0 ] = Latt.dims[ ND - 1 ] ;
    fwrite( in , sizeof( uint32_t ) , 1 , outfile ) ;
    fwrite( tcorr[ munu ] , sizeof( double ) , Latt.dims[ ND - 1 ] , 
	    outfile ) ; 
  }
  fclose( outfile ) ;

  return ;
}
