/**
   @file spinor_ops.c
   @brief gamma multiplies, spinor ops and adjoints
 */
#include "common.h"

#include "matrix_ops.h" // NC*NC multiplies, daggers and traces
#include "spinor_ops.h" // so I can alphabetise the functions

// conjugate transpose of dirac indices
void
adjoint_spinor( struct spinor *adj ,
		const struct spinor S )
{
  int d1 , d2 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      dagger_gauge( (double complex*)adj -> D[d2][d1].C ,
		    (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// trace of the product of two spinors
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B )
{
  int d1 , d2 ;
  register double complex sum = 0.0 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      sum += colortrace_prod( (double complex*)A.D[d1][d2].C ,
			      (double complex*)B.D[d2][d1].C ) ;
    }
  }
  return sum ;
}

// computes G5 ( adj( S ) ) G5
void
full_adj( struct spinor *adj ,
	  const struct spinor S ,
	  const struct gamma G5 )
{
  struct spinor tmp = S ;
  gamma_mul_l( &tmp , G5 ) ;     // left multiply by gamma_5
  gamma_mul_r( &tmp , G5 ) ;     // right multiply by gamma_5
  adjoint_spinor( adj , tmp ) ;  // daggers a spinor
  return ;
}

// atomic left multiply by a gamma matrix
void
gamma_mul_l( struct spinor *res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  int i , j , c1 , c2 ;
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col = GAMMA.ig[i] ;
    double complex fac = 1.0 ;
    // switch for the phases
    switch( GAMMA.g[i] ) {
    case 0 : fac =  1 ; break ;
    case 1 : fac =  I ; break ;
    case 2 : fac = -1 ; break ;
    case 3 : fac = -I ; break ;
    }
    // multiply out
    for( j = 0 ; j < NS ; j++ ) {
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)res -> D[col][j].C ) ;
    } 
    //
  }
  //
  *res = tmp ;
  return ;
}

// multiply a spinor on the right with a gamma matrix
void
gamma_mul_r( struct spinor *res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  int i , j , c1 , c2 ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const int col = GAMMA.ig[j] ;
    double complex fac = 1.0 ;
    switch( GAMMA.g[ col ] ) {
    case 0 : fac =  1 ; break ;
    case 1 : fac =  I ; break ;
    case 2 : fac = -1 ; break ;
    case 3 : fac = -I ; break ;
    }
    // and copy it in
    for( i = 0 ; i < NS ; i++ ) {
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)res -> D[i][col].C ) ;
    } 
    //
  }
  *res = tmp ;
  return ;
}

// multiply by a link :: res = link * S
void
gauge_spinor( struct spinor *res ,  
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
gaugedag_spinor( struct spinor *res ,
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

// multiply by a link :: res = S * link
void
spinor_gauge( struct spinor *res ,  
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

// right multiply link by a daggered spinor res = link * S^{\dagger}
void
gauge_spinordag( struct spinor *res ,
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

// multiply by a daggered link res = S^{\dagger} link
void
spinordag_gauge( struct spinor *res ,
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
spinor_gaugedag( struct spinor *res ,
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
