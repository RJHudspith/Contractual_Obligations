/**
   @file spinor_ops.c
   @brief gauge*spinor, spinor*gauge and daggered variants
 */
#include "common.h"

#include "matrix_ops.h" // NC*NC multiplies, daggers, sums and traces
#include "spinor_ops.h" // so that I can alphabetise

// atomically add two spinors
void
add_spinors( double complex *__restrict SUM ,
	     const double complex *__restrict S )
{
  int i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    SUM[ i ] += S[ i ] ;
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

// sums a propagator over a timeslice
void
sumprop( struct spinor *__restrict SUM ,
	 struct spinor *__restrict S )
{
  zero_spinor( (double complex*)SUM -> D ) ;
  int d1d2 ;
#pragma omp parallel for private(d1d2)
  for( d1d2 = 0 ; d1d2 < NSNS ; d1d2++ ) {
    const int d1 = d1d2 / NS ;
    const int d2 = d1d2 % NS ;
    // accumulate into a matrix
    double complex sum[ NCNC ] = {} ;
    int i ;
    for( i = 0 ; i < VOL3 ; i++ ) {
      add_mat( sum , (double complex*)S[i].D[ d1 ][ d2 ].C ) ;
    }
    colormatrix_equiv( (double complex*)SUM -> D[ d1 ][ d2 ].C , sum ) ;
  }
  return ;
}

// zeros a spinor
void
zero_spinor( double complex *__restrict S )
{
  int i ;
  for( i = 0 ; i < ( NSNS * NCNC ) ; i++ ) {
    S[ i ] = 0.0 ;
  }
  return ;
}
