/**
   @file unit.c
   @brief very simple unit testing framework
 */

#include "common.h"

#include "gammas.h"
#include "minunit.h"     // minimal unit testing framework
#include "spinor_ops.h"  // spinor operations

static struct gamma *GAMMAS = NULL ; // gamma matrix technology

struct spinor A ; // spinor storage

struct latt_info Latt ; // lattice information

// counters
int tests_run = 0 ;
int tests_fail = 0 ;

#define FTOL ( NC * 1.E-14 ) 

static char 
*equate_minus_test( void )
{
  struct spinor B ;
  equate_spinor_minus( &B , (const void*)A.D ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor ops equate minus broken " , 
	       !(*a != -*b) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

static char
*flipsign_test( void )
{
  struct spinor B = A ;
  flipsign_spinor( &B ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  int i ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : spinor ops equate minus broken " , 
	       !(*a != -*b) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

static char 
*gammas_test( void )
{
  mu_assert("[UNIT] error : gamma matrix setup failure ", 
	    !( make_gammas( GAMMAS , NREL ) == FAILURE ) ) ;
  return NULL ;
}

// multiplies spinor by identity matrix and compares
static char
*gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  gauge_spinor( &B , Id , A ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  for( i = 0 ; i < 2*NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : gauge_spinor broken " , 
	       !( *a != *b ) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

// multiply by (-I)^{\dagger}
static char
*gaugedag_test( void )
{
  double complex mI[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      mI[ j + i * NC ] = ( j != i ) ? 0.0 : -I ;
    }
  }
  struct spinor B ;
  gaugedag_spinor( &B , mI , A ) ;
  const double complex *a = (const double complex*)A.D ;
  const double complex *b = (const double complex*)B.D ;
  for( i = 0 ; i < NSNS*NCNC ; i++ ) {
    mu_assert( "[UNIT] error : gaugedag_spinor broken " , 
	       !( fabs( creal( b[i] ) + cimag( a[i] ) ) > FTOL ||
		 fabs( cimag( b[i] ) - creal( a[i] ) ) > FTOL ) ) ;
  }
  return NULL ;
}

// Id * ( A )^{\dagger} where the dagger is JUST OVER COLOR INDICES!
static char
*gauge_spinordag_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  gauge_spinordag( &B , Id , A ) ;
  // B should be the conjugate of the full 12x12 spinor matrix
  int d1 , d2 ;
  for( d1 = 0 ; d1 < NS ; d1++ ) {
    for( d2 = 0 ; d2 < NS ; d2++ ) {
      const double complex *a = (const double complex*)A.D[d1][d2].C ;
      const double complex *b = (const double complex*)B.D[d1][d2].C ;
      int c1 , c2 ;
      for( c1 = 0 ; c1 < NC ; c1++ ) {
	for( c2 = 0 ; c2 < NC ; c2++ ) {			     
	  mu_assert( "[UNIT] error : gauge_spinordag broken " , 
		     !( fabs( creal( b[ c2 + c1*NC ] ) - creal( a[ c1 + NC*c2 ] ) ) > FTOL ||
			fabs( cimag( b[ c2 + c1*NC ] ) + cimag( a[ c1 + NC*c2 ] ) ) > FTOL) ) ;
	}
      }
      //
    }
  }
  return NULL ;
}

// compute B = A * Id, where Id is a color matrix, A & B are spinors
static char
*spinor_gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  spinor_gauge( &B , A , Id ) ;
  const double *a = (const double*)A.D ;
  const double *b = (const double*)B.D ;
  for( i = 0 ; i < 2*NCNC*NSNS ; i++ ) {
    mu_assert( "[UNIT] error : spinor_gauge broken " , 
	       !( *a != *b ) ) ;
    a++ , b++ ;
  }
  return NULL ;
}

// computes B = ( A )^{\dagger} Id, product over color indices
static char
*spinordag_gauge_test( void )
{
  double complex Id[ NCNC ] ;
  int i , j ; 
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      Id[ j + i * NC ] = ( j != i ) ? 0.0 : 1.0 ;
    }
  }
  struct spinor B ;
  spinordag_gauge( &B , A , Id ) ;
  const double complex *a = (const double complex*)A.D ;
  const double complex *b = (const double complex*)B.D ;

  return NULL ;
}

// spinor tests
static char *
spinops_test( void )
{
  // initialise A
  double complex *b = (double complex*)A.D ;
  int i ;
  for( i = 0 ; i < NSNS * NCNC ; i++ ) {
    *b = ( i + I * i ) , b++ ;
  }

  mu_run_test( gammas_test ) ;
  mu_run_test( equate_minus_test ) ;
  mu_run_test( flipsign_test ) ;
  mu_run_test( gauge_test ) ;
  mu_run_test( gaugedag_test ) ;
  mu_run_test( gauge_spinordag_test ) ;
  mu_run_test( spinor_gauge_test ) ;

  return NULL ;
}

int 
main( const int argc , const char *argv[] )
{
  // precompute the gamma basis
  GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // initial gamma setup and test
  char *spinres = spinops_test( ) ;

  if( tests_fail == 0 ) {
    printf( "[UNIT] all %d tests passed \n" ,
	    tests_run ) ;
  } else {
    printf( "%s \n" , spinres ) ;
    printf( "[UNIT] %d out of %d tests failed\n" , 
	    tests_fail , tests_run ) ;
  }

  free( GAMMAS ) ;

  return SUCCESS ;
}
