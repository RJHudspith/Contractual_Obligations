/**
   @file spinor_ops.c
   @brief gamma multiplies, spinor ops and adjoints
 */
#include "common.h"

#include "gammas.h"     // for the gamma_mmuls in meson_contract
#include "matrix_ops.h" // NC*NC multiplies, daggers and traces
#include "spinor_ops.h" // so I can alphabetise the functions

// conjugate transpose of dirac indices
void
adjoint_spinor( struct spinor *__restrict adj ,
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

// convert chiral to non-relativistic
void
chiral_to_nrel( struct spinor *S )
{
  struct spinor P = *S ;
  int c1 , c2 ;
  for( c1 = 0 ; c1 < NC ; c1++ ) {
    for( c2 = 0 ; c2 < NC ; c2++ ) {
      S -> D[0][0].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[1][0].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
      S -> D[2][0].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] + P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[3][0].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] + P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

      S -> D[0][1].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[1][1].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
      S -> D[2][1].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] + P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[3][1].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] + P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );

      S -> D[0][2].C[c1][c2] = 0.5 * ( -P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] + P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[1][2].C[c1][c2] = 0.5 * ( -P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] + P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );
      S -> D[2][2].C[c1][c2] = 0.5 * (  P.D[0][0].C[c1][c2] - P.D[2][0].C[c1][c2] - P.D[0][2].C[c1][c2] + P.D[2][2].C[c1][c2] );
      S -> D[3][2].C[c1][c2] = 0.5 * (  P.D[1][0].C[c1][c2] - P.D[3][0].C[c1][c2] - P.D[1][2].C[c1][c2] + P.D[3][2].C[c1][c2] );

      S -> D[0][3].C[c1][c2] = 0.5 * ( -P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] + P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[1][3].C[c1][c2] = 0.5 * ( -P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] + P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );
      S -> D[2][3].C[c1][c2] = 0.5 * (  P.D[0][1].C[c1][c2] - P.D[2][1].C[c1][c2] - P.D[0][3].C[c1][c2] + P.D[2][3].C[c1][c2] );
      S -> D[3][3].C[c1][c2] = 0.5 * (  P.D[1][1].C[c1][c2] - P.D[3][1].C[c1][c2] - P.D[1][3].C[c1][c2] + P.D[3][3].C[c1][c2] );	
    }
  }
  return ;
}

// computes G5 ( adj( S ) ) G5
void
full_adj( struct spinor *__restrict adj ,
	  const struct spinor S ,
	  const struct gamma G5 )
{
  struct spinor tmp = S ;
  gamma_mul_lr( &tmp , G5 , G5 ) ;     // left multiply by gamma_5
  adjoint_spinor( adj , tmp ) ;  // daggers a spinor
  return ;
}

// atomic left multiply by a gamma matrix
void
gamma_mul_l( struct spinor *__restrict res ,
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
gamma_mul_r( struct spinor *__restrict res ,
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

//
void
gamma_mul_lr( struct spinor *__restrict S , 
	      const struct gamma GLEFT ,
	      const struct gamma GRIGHT )
{
  struct spinor tmp = *S ; // temporary space
  double complex fac = 1.0 ;
  int i , j ;
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col1 = GLEFT.ig[ i ] ;
    for( j = 0 ; j < NS ; j++ ) {
      const int col2 = GRIGHT.ig[ j ] ;
      // switch for the phases
      switch( ( GLEFT.g[ i ] + GRIGHT.g[ col2 ] ) & 3 ) {
      case 0 : fac =  1 ; break ;
      case 1 : fac =  I ; break ;
      case 2 : fac = -1 ; break ;
      case 3 : fac = -I ; break ;
      }
      // multiply out
      constant_mul_gauge( (double complex*)tmp.D[i][j].C , fac ,
			  (const double complex*)S -> D[col1][col2].C ) ;
    }
  }
  //
  *S = tmp ;
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

// meson contraction code
double complex
meson_contract( const struct gamma GSNK ,		
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 )
{
  struct gamma G1 , G2 ;
  gamma_mmul( &G2 , G5 , GSRC ) ;
  gamma_mmul( &G1 , GSNK , G5 ) ;

  double complex sum = 0.0 ;
  register double sumr , sumi ;

  int i , j , c1 , c2 ;
  for( j = 0 ; j < NS ; j++ ) {

    const int col2 = G2.ig[ j ] ;

    // loop columns
    for( i = 0 ; i < NS ; i++ ) {
	
      const int col1 = G1.ig[ i ] ;
	
      // sums in double to avoid complex multiply
      sumr = sumi = 0.0 ;
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  sumr += creal( bwd.D[col2][col1].C[c2][c1] ) * creal( fwd.D[j][i].C[c2][c1] ) 
	    + cimag( bwd.D[col2][col1].C[c2][c1] ) * cimag( fwd.D[j][i].C[c2][c1] ) ;
	  sumi += creal( bwd.D[col2][col1].C[c2][c1] ) * cimag( fwd.D[j][i].C[c2][c1] ) -
	    cimag( bwd.D[col2][col1].C[c2][c1] ) * creal( fwd.D[j][i].C[c2][c1] ) ;
	}
      }

      // switch for the phases
      switch( ( G1.g[ i ] + G2.g[ col2 ] ) & 3 ) {
      case 0 : sum +=  sumr + I * sumi ; break ;
      case 1 : sum += -sumi + I * sumr ; break ;
      case 2 : sum += -sumr - I * sumi ; break ;
      case 3 : sum +=  sumi - I * sumr ; break ;
      }
      // and we are done
    }
  }
  return sum ;
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


