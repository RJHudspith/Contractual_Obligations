/**
   @file contractions.c
   @brief contraction codes and brute force gamma multiplies
 */

#include "common.h"

#include "contractions.h" // so we can alphabetise
#include "matrix_ops.h"   // colortrace_prod

// conjugate transpose of dirac indices
void
adjoint_spinor( struct spinor *__restrict adj ,
		const struct spinor S )
{
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      dagger_gauge( (double complex*)adj -> D[d2][d1].C ,
		    (const double complex*)S.D[d1][d2].C ) ;
    }
  }
  return ;
}

// returns spin-color trace : Tr[ A B ]
double complex
bilinear_trace( const struct spinor A ,
		const struct spinor B )
{
  int d1 , d2 ;
  register double complex sum = 0.0 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      sum += colortrace_prod( (const double complex*)A.D[d1][d2].C ,
			      (const double complex*)B.D[d2][d1].C ) ;
    }
  }
  return sum ;
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
  int i , j ;
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
  int i , j ;
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

// meson contraction code computes Tr[ GSNK ( G5 bwd G5 )^{\dagger} GSRC ( fwd ) ]
double complex
meson_contract( const struct gamma GSNK ,		
		const struct spinor bwd , 
		const struct gamma GSRC ,
		const struct spinor fwd ,
		const struct gamma G5 )
{
  register double gsumr = 0.0 , gsumi = 0.0 ;

  int i , j , c1 , c2 , col2 , col1 ;
  for( j = 0 ; j < NS ; j++ ) {
    
    col2 = GSRC.ig[ G5.ig[ j ] ] ; 

    // loop columns
    for( i = 0 ; i < NS ; i++ ) {

      col1 = G5.ig[ GSNK.ig[ i ] ] ;
      
      // sums in double to avoid complex multiply
      register double sumr = 0.0 , sumi = 0.0 ;
      for( c2 = 0 ; c2 < NC ; c2++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  sumr += creal( bwd.D[col2][col1].C[c2][c1] ) * creal( fwd.D[j][i].C[c2][c1] )
	    + cimag( bwd.D[col2][col1].C[c2][c1] ) * cimag( fwd.D[j][i].C[c2][c1] ) ;
	  sumi += creal( bwd.D[col2][col1].C[c2][c1] ) * cimag( fwd.D[j][i].C[c2][c1] ) -
	    cimag( bwd.D[col2][col1].C[c2][c1] ) * creal( fwd.D[j][i].C[c2][c1] ) ;
	}
      }

      // switch for the phases
      switch( ( GSNK.g[ i ] + G5.g[ col1 ] + G5.g[ col2 ] + GSRC.g[ col2 ] ) & 3 ) {
      case 0 : gsumr +=  sumr ; gsumi +=  sumi ; break ;
      case 1 : gsumr += -sumi ; gsumi +=  sumr ; break ;
      case 2 : gsumr += -sumr ; gsumi += -sumi ; break ;
      case 3 : gsumr +=  sumi ; gsumi += -sumr ; break ;
      }
      // and we are done
    }
  }
  return gsumr + I * gsumi ;
}

// meson contraction code computes Tr[ GSNK ( bwd ) GSRC ( fwd ) ]
double complex
simple_meson_contract( const struct gamma GSNK ,		
		       const struct spinor bwd , 
		       const struct gamma GSRC ,
		       const struct spinor fwd )
{
  register double gsumr = 0.0 , gsumi = 0.0 ;

  int i , j , c1 , c2 ;
  for( i = 0 ; i < NS ; i++ ) {

    const int col1 = GSNK.ig[ i ] ;

    // loop columns
    for( j = 0 ; j < NS ; j++ ) {

      const int col2 = GSRC.ig[ j ] ;

      register double sumr = 0.0 , sumi = 0.0 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  sumr += creal( bwd.D[col1][col2].C[c1][c2] ) * creal( fwd.D[j][i].C[c2][c1] ) 
	    - cimag( bwd.D[col1][col2].C[c1][c2] ) * cimag( fwd.D[j][i].C[c2][c1] ) ;
	  sumi += creal( bwd.D[col1][col2].C[c1][c2] ) * cimag( fwd.D[j][i].C[c2][c1] ) 
	    + cimag( bwd.D[col1][col2].C[c1][c2] ) * creal( fwd.D[j][i].C[c2][c1] ) ;
	}
      }
      // switch for the phases
      switch( ( GSNK.g[ i ] + GSRC.g[ col2 ] ) & 3 ) {
      case 0 : gsumr +=  sumr ; gsumi +=  sumi ; break ;
      case 1 : gsumr += -sumi ; gsumi +=  sumr ; break ;
      case 2 : gsumr += -sumr ; gsumi += -sumi ; break ;
      case 3 : gsumr +=  sumi ; gsumi += -sumr ; break ;
      }
      //
    }
  }
  return gsumr + I * gsumi ;
}
