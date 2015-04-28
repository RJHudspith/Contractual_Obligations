/**
   @file spinor_ops.c
   @brief gauge*spinor, spinor*gauge and daggered variants
 */
#include "common.h"

#include "matrix_ops.h" // NC*NC multiplies, daggers, sums and traces
#include "spinor_ops.h" // so that I can alphabetise

#ifndef HAVE_EMMINTRIN_H

// atomically add two spinors
static void
add_spinors( double complex *SUM ,
	     const double complex *S ) 
{
  int i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    *SUM += *S ; SUM++ ; S++ ;
  }
  return ;
}

// inline zero a spinor
static void
zero_spinor( double complex *S ) 
{
  int i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    *S = 0.0 ; S++ ;
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
  int i ;
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
  int i ;
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
  int d1 , d2 ;
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
  int d1 , d2 ;
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
  int d1 , d2 ;
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

// multiply by a link :: res = S * link
void
spinor_gauge( struct spinor *__restrict res ,  
	      const struct spinor S ,
	      const double complex link[ NCNC ] )
{
  int d1 , d2 ;
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
  int d1 , d2 ;
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
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      multab_dag( (double complex*)res -> D[d1][d2].C , 
		  (const double complex*)S.D[d1][d2].C ,
		  link ) ;
    }
  }
  return ;
}

// sums a propagator over some volume into spinor "SUM"
void
sumprop( void *SUM ,
	 const void *S )
{
  double complex *sum = (double complex*)SUM ;
  zero_spinor( sum ) ;
  const double complex *s = (const double complex*)S ;
  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    add_spinors( sum , s ) ; s += NSNS*NCNC ;
  }
  return ;
}

// multiplies two spinors A = B * A
void
spinmul_atomic_left( struct spinor *A ,
		     const struct spinor B )
{
  struct spinor tmp = *A ;
  int d1 , d2 , d3 ;
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

// zero a spinor
void
spinor_zero( void *S )
{
  struct spinor *s = (struct spinor*)S ;
  int i ;
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

#endif
