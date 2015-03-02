/**
   @file matrix_ops.c
   @brief simple matrix multiplies
 */

#include "common.h"

#include "gammas.h"

// simple NxN square matrix multiplication a = b.c
void 
multab( double complex a[ NCNC ] , 
	const double complex b[ NCNC ] , 
	const double complex c[ NCNC ] )
{
#if NC==3
  a[0] = b[0] * c[0] + b[1] * c[3] + b[2] * c[6] ;	\
  a[1] = b[0] * c[1] + b[1] * c[4] + b[2] * c[7] ;	\
  a[2] = b[0] * c[2] + b[1] * c[5] + b[2] * c[8] ;	\
  a[3] = b[3] * c[0] + b[4] * c[3] + b[5] * c[6] ;	\
  a[4] = b[3] * c[1] + b[4] * c[4] + b[5] * c[7] ;	\
  a[5] = b[3] * c[2] + b[4] * c[5] + b[5] * c[8] ;	\
  a[6] = b[6] * c[0] + b[7] * c[3] + b[8] * c[6] ;	\
  a[7] = b[6] * c[1] + b[7] * c[4] + b[8] * c[7] ;	\
  a[8] = b[6] * c[2] + b[7] * c[5] + b[8] * c[8] ;	
#elif NC==2
  a[0] = b[0] * c[0] + b[1] * c[2] ;		\
  a[1] = b[0] * c[1] + b[1] * c[3] ;		\
  a[2] = b[2] * c[0] + b[3] * c[2] ;		\
  a[3] = b[2] * c[1] + b[3] * c[3] ;		
#else
  // slow and stupid version
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++  ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ j + m*NC ] ) ; IMC = cimag( c[ j + m*NC ] ) ;
	sum += REB * REC - IMB * IMC + I * ( REB * IMC + IMB * REC ) ; 
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}

// 3x3 mult a = ( b^{\dagger} ).c 
void 
multabdag( double complex a[ NCNC ] , 
	   const double complex b[ NCNC ] , 
	   const double complex c[ NCNC ] )
{
#if NC==3
  a[0] = conj( b[0] ) * c[0] + conj( b[3] ) * c[3] + conj( b[6] ) * c[6] ; \
  a[1] = conj( b[0] ) * c[1] + conj( b[3] ) * c[4] + conj( b[6] ) * c[7] ; \
  a[2] = conj( b[0] ) * c[2] + conj( b[3] ) * c[5] + conj( b[6] ) * c[8] ; \
  a[3] = conj( b[1] ) * c[0] + conj( b[4] ) * c[3] + conj( b[7] ) * c[6] ; \
  a[4] = conj( b[1] ) * c[1] + conj( b[4] ) * c[4] + conj( b[7] ) * c[7] ; \
  a[5] = conj( b[1] ) * c[2] + conj( b[4] ) * c[5] + conj( b[7] ) * c[8] ; \
  a[6] = conj( b[2] ) * c[0] + conj( b[5] ) * c[3] + conj( b[8] ) * c[6] ; \
  a[7] = conj( b[2] ) * c[1] + conj( b[5] ) * c[4] + conj( b[8] ) * c[7] ; \
  a[8] = conj( b[2] ) * c[2] + conj( b[5] ) * c[5] + conj( b[8] ) * c[8] ;
#elif NC==2
  a[0] = conj( b[0] ) * c[0] + conj( b[2] ) * c[2] ;	\
  a[1] = conj( b[0] ) * c[1] + conj( b[2] ) * c[3] ;	\
  a[2] = conj( b[1] ) * c[0] + conj( b[3] ) * c[2] ;	\
  a[3] = conj( b[1] ) * c[1] + conj( b[3] ) * c[3] ;
#else
  int i , j , m ;
  register double complex sum ;
  register double REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ i + NC*m ] ) ; IMB = cimag( b[ i + NC*m ] ) ;
	REC = creal( c[ j + NC*m ] ) ; IMC = cimag( c[ j + NC*m ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REB * IMC - IMB * REC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
  return ;
}

// a = b * c^{\dagger}
void 
multab_dag( double complex a[ NCNC ] , 
	    const double complex b[ NCNC ] , 
	    const double complex c[ NCNC ] )
{
#if NC==3
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; \
  a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; \
  a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; \
  a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; \
  a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; \
  a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ; \
  a[6] = b[6] * conj( c[0] ) + b[7] * conj( c[1] ) + b[8] * conj( c[2] ) ; \
  a[7] = b[6] * conj( c[3] ) + b[7] * conj( c[4] ) + b[8] * conj( c[5] ) ; \
  a[8] = b[6] * conj( c[6] ) + b[7] * conj( c[7] ) + b[8] * conj( c[8] ) ; 
#elif NC==2
  a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) ;	\
  a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) ;	\
  a[2] = b[2] * conj( c[0] ) + b[3] * conj( c[1] ) ;	\
  a[3] = b[2] * conj( c[2] ) + b[3] * conj( c[3] ) ;
#else // instead of inlining we have a function call
  int i , j , m ;
  register GLU_complex sum ;
  register GLU_real REB , IMB , REC , IMC ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      sum = 0. ;
      for( m = 0 ; m < NC ; m++ ) {
	REB = creal( b[ m + NC*i ] ) ; IMB = cimag( b[ m + NC*i ] ) ;
	REC = creal( c[ m + NC*j ] ) ; IMC = cimag( c[ m + NC*j ] ) ;
	sum += REB * REC + IMB * IMC + I * ( REC * IMB - REB * IMC ) ;
      }
      a[ j + NC*i ] = sum ;
    }
  }
#endif
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

// is just Tr( a * b )
static double complex
colortrace_prod( double complex *a , double complex *b )
{
  return a[0] * b[0] + a[1] * b[3] + a[2] * b[6] +	\
    a[3] * b[1] + a[4] * b[4] + a[5] * b[7] +		\
    a[6] * b[2] + a[7] * b[5] + a[8] * b[8] ;
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

static void
write_full_spinor( const struct spinor S )
{
  int d1 , d2 , c1 , c2 ;
  printf( "\n" ) ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( c1 = 0 ; c1 < NC ; c1++ ) {

      for( d2 = 0 ; d2 < NS ; d2++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  printf( " %1.2f %1.2f " , 
		  creal( S.D[d1][d2].C[c1][c2] ) ,
		  cimag( S.D[d1][d2].C[c1][c2] ) ) ;
	}
      }
      printf( "\n" ) ;
    }
  }

  return ;
}

// left multiply s.t. res = GAMMA * src
void
gamma_mul_l( struct spinor *res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  int i , j , c1 , c2 ;
  // loop columns
  for( i = 0 ; i < NS ; i++ ) {
    const int col = GAMMA.ig[i] ;
    // switch for the phases
    switch( GAMMA.g[i] ) {
    case 0 : 
      for( j = 0 ; j < NS ; j++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = res -> D[col][j].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    case 1 : 
      for( j = 0 ; j < NS ; j++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = I * res -> D[col][j].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    case 2 : 
      for( j = 0 ; j < NS ; j++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = -res -> D[col][j].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    case 3 : 
      for( j = 0 ; j < NS ; j++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = -I * res -> D[col][j].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    }
  }
  //
  *res = tmp ;
  return ;
}

// atomic multiply of a gamma matrix
void
gamma_mul_r( struct spinor *res ,
	     const struct gamma GAMMA )
{
  struct spinor tmp = *res ; // temporary space
  int i , j , c1 , c2 ;
  // loop columns of src
  for( j = 0 ; j < NS ; j++ ) {
    const int col = GAMMA.ig[j] ;
    switch( GAMMA.g[ col ] ) {
    case 0 :
      for( i = 0 ; i < NS ; i++ ) { 
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = res -> D[i][col].C[c1][c2] ; 
	  }
	}
      } break ;
    case 1 : 
      for( i = 0 ; i < NS ; i++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = I * res -> D[i][col].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    case 2 : 
      for( i = 0 ; i < NS ; i++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = -res -> D[i][col].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    case 3 : 
      for( i = 0 ; i < NS ; i++ ) {
	for( c1 = 0 ; c1 < NC ; c1++ ) {
	  for( c2 = 0 ; c2 < NC ; c2++ ) {
	    tmp.D[i][j].C[c1][c2] = -I * res -> D[i][col].C[c1][c2] ; 
	  }
	}
      } 
      break ;
    }
  }
  *res = tmp ;
  return ;
}

// conjugate transpose of dirac indices
void
adjoint_spinor( struct spinor *adj ,
		const struct spinor S )
{
  int d1 , d2 , c1 , c2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      //
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {
	  adj -> D[d2][d1].C[c2][c1] = conj( S.D[d1][d2].C[c1][c2] ) ;
	}
      }
      //
    }
  }
  return ;
}

// take the adjoint of a gamma matrix
void
adjoint_gamma( struct gamma *adj ,
	       const struct gamma G )
{
  int i ;
  for( i = 0 ; i < NS ; i++ ) {
    const int col = G.ig[i] ;
    adj -> ig[col] = G.ig[ col ] ;
    switch( G.g[i] ) {
    case 0 : adj -> g[col] = 0 ; break ;
    case 1 : adj -> g[col] = 3 ; break ;
    case 2 : adj -> g[col] = 2 ; break ;
    case 3 : adj -> g[col] = 1 ; break ;
    }
  }
  return ;
}

// computes G5 ( adj( S ) ) G5
void
full_adj( struct spinor *adj ,
	  const struct spinor S ,
	  const struct gamma G5 )
{
  struct spinor tmp = S ;
  gamma_mul_l( &tmp , G5 ) ;    // left multiply by gamma_5
  gamma_mul_r( &tmp , G5 ) ;    // right multiply by gamma_5
  adjoint_spinor( adj , tmp ) ;  // daggers a spinor
  return ;
}
