/**
   @file halfspinor_ops.c
   @brief perform arithmetic operations on halfspinors
 */
#include "common.h"

#include "matrix_ops.h"

// set s-matrix to zero
void
zero_halfspinor( struct halfspinor *S )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    zero_colormatrix( S[i].D[0] ) ;
    zero_colormatrix( S[i].D[1] ) ;
    zero_colormatrix( S[i].D[2] ) ;
    zero_colormatrix( S[i].D[3] ) ;
  }
  return ;
}

void
halfspinor_Saxpy( struct halfspinor *H ,
		  const struct halfspinor *S ,
		  const double fac )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_Saxpy( H[i].D[0] , S[i].D[0] , fac ) ;
    colormatrix_Saxpy( H[i].D[1] , S[i].D[1] , fac ) ;
    colormatrix_Saxpy( H[i].D[2] , S[i].D[2] , fac ) ;
    colormatrix_Saxpy( H[i].D[3] , S[i].D[3] , fac ) ;
  }
  return ;
}

// version multiplied by I
void
halfspinor_iSaxpy( struct halfspinor *H ,
		   const struct halfspinor *S ,
		   const double fac )
{
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    colormatrix_iSaxpy( H[i].D[0] , S[i].D[0] , fac ) ;
    colormatrix_iSaxpy( H[i].D[1] , S[i].D[1] , fac ) ;
    colormatrix_iSaxpy( H[i].D[2] , S[i].D[2] , fac ) ;
    colormatrix_iSaxpy( H[i].D[3] , S[i].D[3] , fac ) ;
  }
  return ;
}

void
halfspinor_sigma_Saxpy( struct halfspinor *H ,
			const struct halfspinor *S ,
			const size_t sigma_map[ NS ] ,
			const uint8_t imap[ NS ] )
{
  size_t i , d ;
  for( d = 0 ; d < NS ; d++ ) {
    switch( imap[d]%NS ) {
    case 0 :
      for( i = 0 ; i < LCU ; i++ ) {
	colormatrix_Saxpy( H[i].D[d] , S[i].D[ sigma_map[d] ] , 1. ) ;
      }
      break ;
    case 1 :
      for( i = 0 ; i < LCU ; i++ ) {
	colormatrix_iSaxpy( H[i].D[d] , S[i].D[ sigma_map[d] ] , 1. ) ;
      }
      break ;
    case 2 :
      for( i = 0 ; i < LCU ; i++ ) {
	colormatrix_Saxpy( H[i].D[d] , S[i].D[ sigma_map[d] ] , -1. ) ;
      }
      break ;
    case 3 :
      for( i = 0 ; i < LCU ; i++ ) {
	colormatrix_iSaxpy( H[i].D[d] , S[i].D[ sigma_map[d] ] , -1. ) ;
      }
      break ;
    }
  }
  return ;
}

void
colormatrix_halfspinor( struct halfspinor *a ,
			const double complex b[ NCNC ] ,
			struct halfspinor c )
{
  multab( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multab( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multab( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multab( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
}

void
halfspinor_colormatrix( struct halfspinor *a ,
			const struct halfspinor b ,
			const double complex c[ NCNC ] )
{
  multab( (void*)a->D[0] , (void*)b.D[0] , (void*)c ) ;
  multab( (void*)a->D[1] , (void*)b.D[1] , (void*)c ) ;
  multab( (void*)a->D[2] , (void*)b.D[2] , (void*)c ) ;
  multab( (void*)a->D[3] , (void*)b.D[3] , (void*)c ) ;
}

void
halfspinor_multiply( struct halfspinor *a ,
		     const struct halfspinor b ,
		     const struct halfspinor c )
{
  double complex temp[ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  // first element
  multab( (void*)a -> D[0] , (void*)b.D[0] , (void*)c.D[0] ) ;
  multab( (void*)temp , (void*)b.D[1] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[0] , (void*)temp ) ;
  // second
  multab( (void*)a -> D[1] , (void*)b.D[0] , (void*)c.D[1] ) ;
  multab( (void*)temp , (void*)b.D[1] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[1] , (void*)temp ) ;
  // third
  multab( (void*)a -> D[2] , (void*)b.D[2] , (void*)c.D[0] ) ;
  multab( (void*)temp , (void*)b.D[3] , (void*)c.D[2] ) ;
  add_mat( (void*)a -> D[2] , (void*)temp ) ;
  // fourth
  multab( (void*)a -> D[3] , (void*)b.D[2] , (void*)c.D[1] ) ;
  multab( (void*)temp , (void*)b.D[3] , (void*)c.D[3] ) ;
  add_mat( (void*)a -> D[3] , (void*)temp ) ;
  return ;
}
