/**
   @file halfspinor_ops.c
   @brief perform arithmetic operations on halfspinors
 */
#include "common.h"

#include "matrix_ops.h"
#include "mmul.h"

#ifndef HAVE_IMMINTRIN_H

// atomically add halfspinors a += b
void
add_halfspinor( struct halfspinor *a ,
		const struct halfspinor b )
{
  add_mat( (void*)a -> D[0] , (void*)b.D[0] ) ;
  add_mat( (void*)a -> D[1] , (void*)b.D[1] ) ;
  add_mat( (void*)a -> D[2] , (void*)b.D[2] ) ;
  add_mat( (void*)a -> D[3] , (void*)b.D[3] ) ;
  return ;
}

void
Fmunu_halfspinor( struct halfspinor *a ,
		  const double complex *b ,
		  const struct halfspinor c )
{
  multab( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multab( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multab( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multab( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
  return ;
}

void
colormatrix_halfspinor( double complex *a ,
			const double complex *b ,
			const double complex *c )
{
  multab( (void*)a , (void*)b , (void*)c ) ; a+=NCNC ; c+=NCNC ;
  multab( (void*)a , (void*)b , (void*)c ) ; a+=NCNC ; c+=NCNC ;
  multab( (void*)a , (void*)b , (void*)c ) ; a+=NCNC ; c+=NCNC ;
  multab( (void*)a , (void*)b , (void*)c ) ;
  return ;
}

void
colormatrixdag_halfspinor( struct halfspinor *a ,
			   const double complex b[ NCNC ] ,
			   const struct halfspinor c )
{
  multabdag( (void*)a->D[0] , (void*)b , (void*)c.D[0] ) ;
  multabdag( (void*)a->D[1] , (void*)b , (void*)c.D[1] ) ;
  multabdag( (void*)a->D[2] , (void*)b , (void*)c.D[2] ) ;
  multabdag( (void*)a->D[3] , (void*)b , (void*)c.D[3] ) ;
  return ;
}

void
halfspinor_Saxpy( struct halfspinor *H ,
		  const struct halfspinor S ,
		  const double fac )
{
  colormatrix_Saxpy( H -> D[0] , S.D[0] , fac ) ;
  colormatrix_Saxpy( H -> D[1] , S.D[1] , fac ) ;
  colormatrix_Saxpy( H -> D[2] , S.D[2] , fac ) ;
  colormatrix_Saxpy( H -> D[3] , S.D[3] , fac ) ;
  return ;
}

// version multiplied by I
void
halfspinor_iSaxpy( struct halfspinor *H ,
		   const struct halfspinor S ,
		   const double fac )
{
  colormatrix_iSaxpy( H -> D[0] , S.D[0] , fac ) ;
  colormatrix_iSaxpy( H -> D[1] , S.D[1] , fac ) ;
  colormatrix_iSaxpy( H -> D[2] , S.D[2] , fac ) ;
  colormatrix_iSaxpy( H -> D[3] , S.D[3] , fac ) ;
  return ;
}

// does H^i = H^i + \sigma^ij S^j
// where sigma_map is an index map for the sigma product
// and imap is the map in Z_4 of the elements.
// For example :
// sigma_map[ NS ] = { 2 , 3 , 0 , 1 }
// as sigma_x is { 0 , 1 , 1 , 0 }
// and imap would be { 0 , 0 , 0 , 0 } as this maps to 1
void
halfspinor_sigma_Saxpy( struct halfspinor *H ,
			const struct halfspinor S ,
			const uint8_t sigma_map[ NS ] ,
			const uint8_t imap[ NS ] )
{
  size_t d ;
  for( d = 0 ; d < NS ; d++ ) {
    switch( imap[d]%NS ) {
    case 0 :
      colormatrix_Saxpy( H -> D[d] , S.D[ sigma_map[d] ] , 1. ) ;
      break ;
    case 1 :
      colormatrix_iSaxpy( H -> D[d] , S.D[ sigma_map[d] ] , 1. ) ;
      break ;
    case 2 :
      colormatrix_Saxpy( H -> D[d] , S.D[ sigma_map[d] ] , -1. ) ;
      break ;
    case 3 :
      colormatrix_iSaxpy( H -> D[d] , S.D[ sigma_map[d] ] , -1. ) ;
      break ;
    }
  }
  return ;
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

// set s-matrix to zero
void
zero_halfspinor( struct halfspinor *S )
{
  zero_colormatrix( S -> D[0] ) ;
  zero_colormatrix( S -> D[1] ) ;
  zero_colormatrix( S -> D[2] ) ;
  zero_colormatrix( S -> D[3] ) ;
  return ;
}

#endif
