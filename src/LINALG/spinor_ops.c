/**
   @file spinor_ops.c
   @brief gauge*spinor, spinor*gauge and daggered variants
 */
#include "common.h"

#include "mmul.h"
#include "matrix_ops.h" // NC*NC multiplies, daggers, sums and traces
#include "spinor_ops.h" // so that I can alphabetise

#ifndef HAVE_EMMINTRIN_H

// sum over a spinor
static void
sum_spinor( double complex *SUM ,
	    const double complex *S ) 
{
  size_t i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    *SUM += *S ; SUM++ ; S++ ;
  }
  return ;
}


// inline zero a spinor
static void
zero_spinor( double complex *S ) 
{
  size_t i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    *S = 0.0 ; S++ ;
  }
  return ;
}

// atomically add spinors
void
add_spinors( struct spinor *A ,
	     const struct spinor B )
{
  double complex *s1 = (double complex*)A -> D ;
  double complex *s2 = (double complex*)B.D ;
  size_t i ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s1 += *s2 ; s1++ ; s2++ ;
  }
  return ;
}

// trace out the color indices, S1 is a NS*NS flattened matrix
void
colortrace_spinor( void *S1 ,
		   const void *S2 )
{
  double complex *s1 = (double complex*)S1 ;
  const struct spinor *s2 = (const struct spinor*)S2 ;
  size_t i , j ;
  for( i = 0 ; i < NS ; i++ ) {
    for( j = 0 ; j < NS ; j++ ) {
      s1[ j + i*NS ] = colortrace( (double complex*)s2->D[i][j].C ) ;
    }
  }
  return ;
}

// equate spinors
void
equate_spinor( void *S1 ,
	       const void *S2 )
{
  // probably be better calling out to memcpy ....
  memcpy( S1 , S2 , sizeof( struct spinor ) ) ;
  return ;
}

// equate one spinor to the minus of another
void
equate_spinor_minus( void *mS ,
		     const void *S )
{
  double complex *m = (double complex*)mS ;
  const double complex *s = (double complex*)S ;
  size_t i ;
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    *m = -*s ; m++ ; s++ ;
  }
  return ;
}

// flip a spinor, call with flipsign_spinor( (double complex*)S.D )
void
flipsign_spinor( void *S ) 
{
  double complex *s = (double complex*)S ;
  size_t i ;
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    *s = -*s ; s++ ;
  }
  return ;
}

// multiply by a link :: res = link * S
void
gauge_spinor( struct spinor *__restrict res ,  
	      const double complex link[ NCNC ] ,
	      const struct spinor S )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      // computes S1 = link * S1
      multab( (double complex*)res -> D[d1][d2].C , 
	      link ,
	      (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// multiply by a daggered link res = link^{\dagger} * S
void
gaugedag_spinor( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      // computes S1 = link^{\dagger} * S1
      multabdag( (double complex*)res -> D[d1][d2].C , 
		 link ,
		 (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// right multiply link by a daggered spinor res = link * S^{\dagger}
void
gauge_spinordag( struct spinor *__restrict res ,
		 const double complex link[ NCNC ] ,
		 const struct spinor S )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      // computes S1 = link^{\dagger} * S1
      multab_dag( (double complex*)res -> D[d1][d2].C , 
		  link ,
		  (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// set a spinor to the identity
void
identity_spinor( struct spinor *__restrict res )
{
  zero_spinor( (double complex*)res -> D ) ;
  size_t d , c ;
  for( d = 0 ; d < NS ; d++ ) {
    for( c = 0 ; c < NC ; c++ ) {
      res -> D[ d ][ d ].C[ c ][ c ] = 1.0 ;
    }
  }
  return ;
}

// atomically subtract two spinors
void
sub_spinors( struct spinor *A ,
	      const struct spinor B )
{
  double complex *s1 = (double complex*)A -> D ;
  double complex *s2 = (double complex*)B.D ;
  size_t i ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s1 -= *s2 ; s1++ ; s2++ ;
  }
  return ;
}

// multiply by a link :: res = S * link
void
spinor_gauge( struct spinor *__restrict res ,  
	      const struct spinor S ,
	      const double complex link[ NCNC ] )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      multab( (double complex*)res -> D[d1][d2].C , 
	      (const double complex*)S.D[d1][d2].C , 
	      link ) ;
    }
  }
  return ;
}

// multiply by a daggered link res = S^{\dagger} link
void
spinordag_gauge( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      multabdag( (double complex*)res -> D[d1][d2].C , 
		 (const double complex*)S.D[d1][d2].C ,
		 link ) ;
    }
  }
  return ;
}

// right multiply by a daggered link res = S * link^{\dagger}
void
spinor_gaugedag( struct spinor *__restrict res ,
		 const struct spinor S ,
		 const double complex link[ NCNC ] )
{
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      multab_dag( (double complex*)res -> D[d1][d2].C , 
		  (const double complex*)S.D[d1][d2].C ,
		  link ) ;
    }
  }
  return ;
}

// computest A = A + S * B
void
spinor_Saxpy( struct spinor *A ,
	      const double S ,
	      const struct spinor B )
{
  double complex *s1 = (double complex*)A -> D ;
  const double complex *s2 = (const double complex*)B.D ;
  size_t i ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    *s1 += S * (*s2) ; s1++ ; s2++ ;
  }
  return ;
}

// multiplies two spinors A = B * A
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B )
{
  struct spinor tmp = *A ;
  size_t d1 , d2 , d3 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      double complex link[ NCNC ] , sum[ NCNC ] ;
      // zero
      for( d3 = 0 ; d3 < NCNC ; d3++ ) { sum[ d3 ] = 0.0 ; }
      // spin-color matrix multiply
      for( d3 = 0 ; d3 < NS ; d3++ ) {
	multab( link ,
		(const double complex*)B.D[d1][d3].C ,
		(const double complex*)tmp.D[d3][d2].C ) ;
	add_mat( sum , link ) ;
      }
      colormatrix_equiv( (double complex*)A -> D[d1][d2].C , sum ) ;
    }
  }
  return ;
}

// multiplies two spinors A = A * B
void
spinmul_atomic_right( struct spinor *A ,
		      const struct spinor B )
{
  struct spinor tmp = *A ; // memcpy a spinor here
  size_t d1 , d2 , d3 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      double complex link[ NCNC ] , sum[ NCNC ] ;
      // zero
      for( d3 = 0 ; d3 < NCNC ; d3++ ) { sum[ d3 ] = 0.0 ; }
      // spin-color matrix multiply
      for( d3 = 0 ; d3 < NS ; d3++ ) {
	multab( link ,
		(const double complex*)tmp.D[d1][d3].C ,
		(const double complex*)B.D[d3][d2].C ) ;
	add_mat( sum , link ) ;
      }
      colormatrix_equiv( (double complex*)A -> D[d1][d2].C , sum ) ;
    }
  }
  return ;
}

// zero a spinor
void
spinor_zero( void *S )
{
  struct spinor *s = (struct spinor*)S ;
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    zero_spinor( (double complex*)s[i].D ) ;
  }
  return ;
}

// zero a spinor
void
spinor_zero_site( void *S )
{
  struct spinor *s = (struct spinor*)S ;
  zero_spinor( (double complex*)s -> D ) ;
}

// spintrace into a color matrix
void
spintrace( void *S ,
	   const void *S2 )
{
  double complex *s = (double complex*)S ;
  const struct spinor *s2 = (const struct spinor*)S2 ;
  size_t i , j , d ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      ///
      register double complex sum = 0.0 ;
      for( d = 0 ; d < NS ; d++ ) {
	sum += s2 -> D[d][d].C[i][j] ;
      }
      s[ j + i * NC ] = sum ;
    }
  }
  return ;
}

// sums a propagator over spatial volume into spinor "SUM"
void
sumprop( void *SUM ,
	 const void *S )
{
  double complex *sum = (double complex*)SUM ;
  zero_spinor( sum ) ;
  const double complex *s = (const double complex*)S ;
  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    sum_spinor( sum , s ) ; s += NSNS*NCNC ;
  }
  return ;
}

// sums propagators over spatial volume into spinor "SUM"
void
sumwalls( struct spinor *SUM ,
	  const struct spinor **S ,
	  const size_t Nprops )
{
  size_t mu ;
#pragma omp parallel for private(mu)
  for( mu = 0 ; mu < Nprops ; mu++ ) {
    sumprop( &SUM[mu] , S[mu] ) ;
  }
  return ;
}

// dirac transpose a spinor, returns S^T on stack
struct spinor
transpose_spinor( const struct spinor S )
{
  struct spinor ST ;
  size_t d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      colormatrix_equiv( (double complex*)ST.D[d2][d1].C ,
			 (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ST ;
}

#endif
