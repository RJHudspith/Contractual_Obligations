/**
   @file PropOps.c
   @brief propagator operations

   Supported ops:

   Prop3 = Prop1 * Prop2
   Prop3 = Prop1 + Prop2
   Prop3 = Prop1 - Prop2

   *multiply will do prop1^{\alpha\beta}_{ab}(x) * prop2^{\beta\kappa}_{bc}(x)
   I am not sure if anyone wants that but we already had the linalg for it

   Prop3 is the outputted propagator
 */
#include "common.h"          // one header to rule them all

#include "geometry.h"        // init_geom()
#include "io.h"              // read_prop
#include "input_reader.h"    // free_props
#include "progress_bar.h"    // progress_bar()
#include "read_propheader.h" // read the propagator file header
#include "spinor_ops.h"      // add_spinors, sub_spinors, spinmul_atomic_left

// lattice information holds dimensions and other stuff
// to be taken from the gauge configuration file OR the input file
struct latt_info Latt ;

struct site *lat = NULL ;

// strings are equal
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// write a propagator header file
static void
write_propheader( FILE *outfile ,
		  const struct propagator prop )
{
  size_t mu ;
  // write out usual header stuff
  fprintf( outfile , "<start_header>\n" ) ;
  fprintf( outfile , "Lattice:" ) ;
  for( mu = 0 ; mu < NS ; mu++ ) {
    fprintf( outfile , " %zu" , Latt.dims[mu] ) ;
  }
  fprintf( outfile , "\n" ) ;
  fprintf( outfile , "SrcPos:" ) ;
  for( mu = 0 ; mu < NS ; mu++ ) {
    fprintf( outfile , " %zu" , prop.origin[mu] ) ;
  }
  fprintf( outfile , "\n" ) ;
  // basis
  fprintf( outfile , "Basis: " ) ;
  switch( prop.basis ) {
  case NREL_FWD : fprintf( outfile , "Nrel_fwd\n" ) ; break ;
  case NREL_BWD : fprintf( outfile , "Nrel_bwd\n" ) ; break ;
  case STATIC : fprintf( outfile , "Nrel_fwd\n" ) ; break ;
  case CHIRAL : fprintf( outfile , "Chiral\n" ) ; break ;
  }
  // source type
  fprintf( outfile , "Source: " ) ;
  switch( prop.source ) {
  case POINT : fprintf( outfile , "Point\n" ) ; break ;
  case WALL : fprintf( outfile , "Wall\n" ) ; break ;
  }
  // Endian-ness
  if( WORDS_BIGENDIAN ) {
    fprintf( outfile , "Endian: Big\n" ) ;
  } else {
    fprintf( outfile , "Endian: Little\n" ) ;
  }
  // Precision will always be single!
  fprintf( outfile , "Precision: Single\n" ) ;
  fprintf( outfile , "<end_header>\n" ) ;
  return ;
}

// enumerate the inputs
enum{ DIMS = 1 , PROP1 = 2 , OP = 3 , PROP2 = 4 , OUTPROP = 5 } ;

// main file
int 
main( const int argc, 
      const char *argv[] )
{
  // usual check for command-line arguments
  if( argc != 6 ) {
    fprintf( stderr , "USAGE :: ./PROPOPS Lx,Ly,Lz,Lt prop1 op prop2 outprop\n"
	     "Where op is either +,- or *\n" ) ;
    return FAILURE ;
  }

  // normal allocations and stuff
  void (*f)( struct spinor *A , const struct spinor B ) ;
  struct spinor *S1 = NULL , *S2 = NULL ;
  struct propagator *prop = NULL ;
  FILE *outfile = NULL ;
  int error_code = SUCCESS ;

  // parse the lattice dimensions
  char *tok = strtok( (char*)argv[ DIMS ] , "," ) , *endptr ;
  size_t mu = 1 ;
  Latt.dims[0] = strtol( tok , &endptr , 10 ) ;
  while( ( tok = strtok( NULL , "," ) ) != NULL ) {
    Latt.dims[mu] = strtol( tok , &endptr , 10 ) ;
    mu++ ;
  }
  init_geom( ) ;

  // parse the operator and create a function pointer
  if( are_equal( argv[ OP ] , "+" ) ) {
    f = add_spinors ;
  } else if( are_equal( argv[ OP ] , "-" ) ) {
    f = sub_spinors ;
  } else if( are_equal( argv[ OP ] , "*" ) ) {
    f = spinmul_atomic_right ;
  } else {
    fprintf( stderr , "[OP] operator %s not recognised\n" , argv[ OP ] ) ;
    error_code = FAILURE ; goto end ;
  }
  
  // allocate the prop struct
  prop = malloc( 2 * sizeof( struct propagator ) ) ;
  if( ( prop[0].file = fopen( argv[ PROP1 ] , "rb" ) ) == NULL ) {
    fprintf( stderr , "[IO] prop file %s not found\n" , argv[ PROP1 ] ) ;
    error_code = FAILURE ; goto end ;
  }
  if( ( prop[1].file = fopen( argv[ PROP2 ] , "rb" ) ) == NULL ) {
    fprintf( stderr , "[IO] prop file %s not found\n" , argv[ PROP2 ] ) ;
    error_code = FAILURE ; goto end ;
  }
  if( read_propheaders( prop , 2 ) == FAILURE ) {
    error_code = FAILURE ; goto end ;
  }

  // check the propagators are equivalent?
  if( prop[0].basis != prop[1].basis ) {
    fprintf( stderr , "[IO] unequal bases for the props\n" ) ;
    error_code = FAILURE ; goto end ;
  }
  if( prop[0].source != prop[1].source ) {
    fprintf( stderr , "[IO] unequal source types for the props\n" ) ;
    error_code = FAILURE ; goto end ;
  }
  if( prop[0].precision != prop[1].precision ) {
    fprintf( stderr , "[IO] unequal precisions for the props\n" ) ;
    error_code = FAILURE ; goto end ;
  }
  for( mu = 0 ; mu < NS ; mu++ ) {
    if( prop[0].origin[mu] != prop[1].origin[mu] ) {
      fprintf( stderr , "[IO] unequal origins for the props\n"
	       "   %zu -> %zu %zu \n" , mu ,
	       prop[0].origin[mu] , prop[1].origin[mu] ) ;
      error_code = FAILURE ; goto end ;
    }
  }
  
  // read S1 and S2
  if( corr_malloc( (void**)&S1 , 16 , LCU * sizeof( struct spinor ) ) != 0 ||
      corr_malloc( (void**)&S2 , 16 , LCU * sizeof( struct spinor ) ) != 0 ) {
    fprintf( stderr , "[PROPOPS] Allocation failure\n" ) ;
    error_code = FAILURE ; goto end ;
  }

  // result will be in S1
  fprintf( stdout , "[IO] writing operation to %s \n" , argv[ OUTPROP ] ) ;
  
  outfile = fopen( argv[ OUTPROP ] , "wb" ) ;
  write_propheader( outfile , prop[0] ) ;

  float complex *s = malloc( LCU*NSNS*NCNC * sizeof( float complex ) ) ;
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {
    
    read_prop( prop[0] , S1 ) ;
    read_prop( prop[1] , S2 ) ;
    size_t site ;
    for( site = 0 ; site < LCU ; site++ ) {
      f( &S1[site] , S2[site] ) ;
    }

    // write this timeslice of prop to file
    const double complex *Sp = (const double complex*)S1 ;
    for( site = 0 ; site < LCU*NSNS*NCNC ; site++ ) {
      s[site] = Sp[site] ;
    }
    fwrite( s , sizeof( float complex ) , LCU*NSNS*NCNC , outfile ) ;

    progress_bar( t , LT ) ;
  }

  free( s ) ;
  
  fclose( outfile ) ;

 end :

  if( S1 != NULL ) {
    free( S1 ) ;
  }
  if( S2 != NULL ) {
    free( S2 ) ;
  }
  if( prop != NULL ) {
    free_props( prop , 2 ) ;
  }
  
  return error_code ;
}
