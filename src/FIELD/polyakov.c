/**
   @file polyakov.c
   @brief polyakov loop calculator
 */
#include "common.h"

#include "geometry.h"
#include "matrix_ops.h"
#include "spinor_ops.h"

// precompute all of the polyakov lines
static double complex **lines = NULL ;

// simple matrix multiplication a = a * b
static void 
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
  size_t i , j , m ;
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
  size_t i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      poly[ j + i * NC ] = ( i == j ) ? 1.0 : 0.0 ;
    }
  }
#endif
  return ;
}

// precompute all of the lines
void
precompute_poly_lines( const struct site *__restrict lat ,
		       const size_t site , 
		       const size_t dir ,
		       const size_t tsrc , // is a timeslice idx
		       const GLU_bool forward ,
		       const size_t length )
{
  // allocate the lines, there will be LVOLUME of these
  corr_malloc( (void**)&lines , ALIGNMENT , LVOLUME ) ;
  size_t i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    corr_malloc( (void**)lines[i] , ALIGNMENT , NCNC ) ;
  }

  // loop sub volume, multiplying the time like links of the previous timeslice
  size_t t , prev = lat[ tsrc*LCU ].back[ ND-1 ] ;
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {
    // compute this timeslice and the previous one
    const size_t idx = lat[ prev ].neighbor[ ND-1 ] ; 
    // can do this in parallel
    size_t j ;
#pragma omp parallel for private(j)
    for( j = 0 ; j < LCU ; j++ ) {
      #ifdef HAVE_IMMINTRIN_H
      multab( (__m128d*)lines[ j + idx ] ,
	      (__m128d*)lines[ j + prev ] ,
	      (__m128d*)lat[ j + idx ].O[ ND-1 ] ) ;
      #else
      multab( lines[ j + idx ] ,
	      lines[ j + prev ] ,
	      lat[ j + idx ].O[ ND-1 ] ) ;
      #endif
    }
    prev = idx ;
  }

  return ;
}

// compute a short polyakov line
static void
small_poly( double complex poly[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t site , 
	    const size_t dir ,
	    const GLU_bool forward ,
	    const size_t length )
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
  // set the spinor to 0
  spinor_zero( S ) ;

  // put it in the top left sector
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    memcpy( S[i].D[0][0].C , lines[i+LCU*t] , NCNC*sizeof(double complex) ) ;
    memcpy( S[i].D[1][1].C , lines[i+LCU*t] , NCNC*sizeof(double complex) ) ;
  }

#if 0
  // compute the direction it is going in
  const double direction = ( forward == GLU_TRUE ) ? 1 : -1 ;

  // set the spinor to 0
  spinor_zero( S ) ;

  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // compute polyakov line
    double complex pline[ NCNC ] ;
    small_poly( pline , lat , i , ND-1 , forward , t ) ;

    size_t j ;
    for( j = 0 ; j < NS ; j++ ) {

      // multiply by 1/2
      double complex res[ NCNC ] ;
      constant_mul_gauge( res , 0.5 , pline ) ;

      // identity matrix
      colormatrix_equiv( (double complex*)S[i].D[j][j].C , res ) ;
      
      // multiply by gamma
      switch( gamma_t.g[j] ) {
      case 0 : constant_mul_gauge( res ,  0.5 * direction , pline ) ; break ;
      case 1 : constant_mul_gauge( res ,  0.5 * I * direction , pline ) ; break ;
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
#endif
}

// computes a polyakov line, looping around the dimension "dir"
double complex
poly( const struct site *__restrict lat , 
      size_t dir )
{
  // if you have stupidly set the dimension to something unreasonable
  // default the direction to ND
  if( dir > ND ) { dir = ND ; }
  double complex sum = 0.0 ;
  size_t i ; 
  //#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LCU ; i++ ) {
    double complex poly[ NCNC ] ;
    // use the correct site for one of the hypercubes ...
    int x[ ND ] ;
    get_mom_2piBZ( x , i , dir ) ;
    const size_t k = gen_site( x ) ;
    small_poly( poly , lat , k , dir , GLU_TRUE , Latt.dims[dir] ) ;
    double complex s ;
    #ifdef HAVE_IMMINTRIN_H
    //_mm_store_pd( (void*)&s , colortrace( (const __m128d*)poly ) ) ;
    #else
    //s = colortrace( poly ) ;
    #endif
    sum = sum + s ;
  }
  return sum ;
}
