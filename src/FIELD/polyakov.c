/**
   @file polyakov.c
   @brief polyakov loop calculator
 */
#include "common.h"

#include "geometry.h"
#include "matrix_ops.h"
#include "spinor_ops.h"

// simple matrix multiplication a = a * b
void 
multab_atomic_right( double complex a[ NCNC ] , 
		     const double complex b[ NCNC ] )
{
#if NC==3
  double complex R0 = a[0] * b[0] + a[1] * b[3] + a[2] * b[6] ;	\
  double complex R1 = a[0] * b[1] + a[1] * b[4] + a[2] * b[7] ;	\
  double complex R2 = a[0] * b[2] + a[1] * b[5] + a[2] * b[8] ;	\
  a[0] = R0 ; a[1] = R1 ; a[2] = R2 ;				\
  R0 = a[3] * b[0] + a[4] * b[3] + a[5] * b[6] ;		\
  R1 = a[3] * b[1] + a[4] * b[4] + a[5] * b[7] ;		\
  R2 = a[3] * b[2] + a[4] * b[5] + a[5] * b[8] ;		\
  a[3] = R0 ; a[4] = R1 ; a[5] = R2 ;				\
  R0 = a[6] * b[0] + a[7] * b[3] + a[8] * b[6] ;		\
  R1 = a[6] * b[1] + a[7] * b[4] + a[8] * b[7] ;		\
  R2 = a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;		\
  a[6] = R0 ; a[7] = R1 ; a[8] = R2 ;				
#elif NC==2
  double complex R0 = a[0] * b[0] + a[1] * b[2] ;		\
  double complex R1 = a[0] * b[1] + a[1] * b[3] ;		\
  a[0] = R0 ; a[1] = R1 ;					\
  R0 = a[2] * b[0] + a[3] * b[2] ;				\
  R1 = a[2] * b[1] + a[3] * b[3] ;				\
  a[2] = R0 ; a[3] = R1 ;				
#else
  // slow and stupid loopy version
  int i , j , m ;
  double complex R[ NC ] ; // this probably needs to be aligned
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j ++ ) {
      R[j] = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) {
        R[ j ] += a[ m + i*NC ] * b[ j + m*NC ] ;
      }	
    }
    // copy back over to a ...
    for( j = 0 ; j < NC ; j ++ ) {
      a[ j + i*NC ] = R[j] ;
    }
  }
#endif
  return ;
}

// identity matrix
static void
identity( double complex poly[ NCNC ] )
{
#if NC == 3
  poly[0] = 1.0 ; poly[1] = 0.0 ; poly[2] = 0.0 ;
  poly[3] = 0.0 ; poly[4] = 1.0 ; poly[5] = 0.0 ;
  poly[6] = 0.0 ; poly[7] = 0.0 ; poly[8] = 1.0 ;
#elif NC == 2
  poly[0] = 1.0 ; poly[1] = 0.0 ;
  poly[2] = 0.0 ; poly[3] = 1.0 ;
#else
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      poly[ j + i * NC ] = ( i == j ) ? 1.0 : 0.0 ;
    }
  }
#endif
  return ;
}

// compute a short polyakov line
static void
small_poly( double complex poly[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const int site , 
	    const int dir ,
	    const GLU_bool forward ,
	    const int length )
{
  // zero length polyakov line is just the identity
  if( length == 0 ) {
    identity( poly ) ;
    return ;
  }
  colormatrix_equiv( poly , lat[ site ].O[dir] ) ;
  // poly loop calculation ...
  size_t next = site , t ;
  for( t = 1 ; t < length ; t++ ) {
    next = ( forward == GLU_TRUE ) ? lat[ next ].neighbor[ dir ] :	\
      lat[ next ].back[ dir ] ;
    multab_atomic_right( poly , lat[ next ].O[ dir ] ) ;
  }
  return ;
}

// projection for static quark ( 1 +/- \gamma_t ) / 2
void
static_quark( struct spinor *S ,
	      const struct gamma gamma_t ,
	      const size_t t ,
	      const GLU_bool forward )
{
  // compute the direction it is going in
  const double direction = ( forward == GLU_TRUE ) ? 1 : -1 ;

  // set the spinor to 0
  spinor_zero( S ) ;

  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // compute polyakov line
    double complex pline[ NCNC ] ;
    small_poly( pline , lat , i , ND-1 , forward , t ) ;

    int j ;
    for( j = 0 ; j < NS ; j++ ) {

      // multiply by 1/2
      double complex res[ NCNC ] ;
      constant_mul_gauge( res , 0.5 , pline ) ;

      // identity matrix
      colormatrix_equiv( (double complex*)S[i].D[j][j].C , res ) ;
      
      // multiply by gamma
      switch( gamma_t.g[j] ) {
      case 0 : constant_mul_gauge( res , 0.5 * direction , pline ) ; break ;
      case 1 : constant_mul_gauge( res , 0.5 * I * direction , pline ) ; break ;
      case 2 : constant_mul_gauge( res , -0.5 * direction , pline ) ; break ;
      case 3 : constant_mul_gauge( res , -0.5 * I * direction , pline ) ; break ;
      }

      // poke it in to the right spot, could be on the diagonal so have to add
      const size_t idx = gamma_t.ig[ j ] ;
      #ifdef HAVE_IMMINTRIN_H
      add_mat( (__m128d*)S[i].D[idx][j].C , (const __m128d*)res ) ;
      #else
      add_mat( (double complex*)S[i].D[idx][j].C , res ) ;
      #endif
    }
  }
  return ;
}

// computes a polyakov line, looping around the dimension "dir"
double complex
poly( const struct site *__restrict lat , 
      int dir )
{
  // if you have stupidly set the dimension to something unreasonable
  // default the direction to ND
  if( dir > ND || dir < 0 ) { dir = ND ; }
  double complex sum = 0.0 ;
  int i ; 
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LCU ; i++ ) {
    double complex poly[ NCNC ] ;
    // use the correct site for one of the hypercubes ...
    int x[ ND ] ;
    get_mom_2piBZ( x , i , dir ) ;
    const int k = gen_site( x ) ;
    small_poly( poly , lat , k , dir , GLU_TRUE , Latt.dims[dir] ) ;
    double complex s ;
    #ifdef HAVE_IMMINTRIN_H
    _mm_store_pd( (void*)&s , colortrace( (const __m128d*)poly ) ) ;
    #else
    s = colortrace( poly ) ;
    #endif
    sum = sum + s ;
  }
  return sum ;
}

