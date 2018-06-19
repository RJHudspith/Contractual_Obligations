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

// enumerate the inputs
enum{ DIMS = 1 , PROP1 = 2 , OP = 3 , PROP2 = 4 , OUTPROP = 5 } ;

// enumerate the propagator operation
typedef enum{ ADD , SUB , MUL } Optype ;

// write a propagator header file
static void
write_propheader( FILE *outfile ,
		  const struct propagator prop1 ,
		  const struct propagator prop2 ,
		  const Optype Op )
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
    fprintf( outfile , " %zu" , prop1.origin[mu] ) ;
  }
  fprintf( outfile , "\n" ) ;
  // basis
  fprintf( outfile , "Basis: " ) ;
  switch( prop1.basis ) {
  case NREL_FWD : fprintf( outfile , "Nrel_fwd\n" ) ; break ;
  case NREL_BWD : fprintf( outfile , "Nrel_bwd\n" ) ; break ;
  case NREL_CORR : fprintf( outfile , "Nrel_fwd\n" ) ; break ;
  case CHIRAL : fprintf( outfile , "Chiral\n" ) ; break ;
  }
  // source type
  fprintf( outfile , "Source: " ) ;
  switch( prop1.source ) {
  case POINT : fprintf( outfile , "Point\n" ) ; break ;
  case WALL : fprintf( outfile , "Wall\n" ) ; break ;
  case Z2_WALL : fprintf( outfile , "Z2_Wall\n" ) ; break ;
  }
  // Endian-ness
  if( WORDS_BIGENDIAN ) {
    fprintf( outfile , "Endian: Big\n" ) ;
  } else {
    fprintf( outfile , "Endian: Little\n" ) ;
  }
  // Precision will always be single!
  fprintf( outfile , "Precision: Single\n" ) ;

  // write out the boundaries if they are different we
  // call them one of the special types, this could change but at
  // the moment it will do
  fprintf( outfile , "Boundaries:" ) ;
  if( ( prop1.bound[ND-1] == PERIODIC &&
	prop2.bound[ND-1] == ANTIPERIODIC ) ||
      ( prop1.bound[ND-1] == ANTIPERIODIC &&
	prop2.bound[ND-1] == PERIODIC ) ) { 
    for( mu = 0 ; mu < NS ; mu++ ) {
      switch( Op ) {
      case ADD : fprintf( outfile , " PplusA" ) ; break ;
      case SUB : fprintf( outfile , " PminusA" ) ; break ;
      case MUL : fprintf( outfile , " PmulA" ) ; break ;
      }
    }
  } else {
    for( mu = 0 ; mu < NS ; mu++ ) {
      switch( prop1.bound[mu] ) {
      case PERIODIC : fprintf( outfile , " periodic" ) ; break ;
      case ANTIPERIODIC : fprintf( outfile , " anti-periodic" ) ; break ;
      case PPLUSA : fprintf( outfile , " PplusA" ) ; break ;
      case PMINUSA : fprintf( outfile , " PminusA" ) ; break ;
      case PMULA : fprintf( outfile , " PmulA" ) ; break ;
      }
    }
  }
  fprintf( outfile , "\n" ) ;

  
  fprintf( outfile , "<end_header>\n" ) ;
  return ;
}

// main file
int 
main( const int argc, 
      const char *argv[] )
{
  // usual check for command-line arguments
  if( argc != 6 ) {
    fprintf( stderr , "USAGE :: ./PROPOPS Lx,Ly,Lz,Lt Peri_Prop op Aperi_Prop outprop\n"
	     "Where op is either +,- or *\n" ) ;
    return FAILURE ;
  }

  // normal allocations and stuff
  struct spinor *S1 = NULL , *S2 = NULL ;
  struct propagator *prop = NULL ;
  FILE *outfile = NULL ;
  int error_code = SUCCESS ;
  Optype Op ; // operator type we will use

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
    Op = ADD ;
  } else if( are_equal( argv[ OP ] , "-" ) ) {
    Op = SUB ;
  } else if( are_equal( argv[ OP ] , "*" ) ) {
    Op = MUL ;
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

  // write the propheader
  outfile = fopen( argv[ OUTPROP ] , "wb" ) ;
  write_propheader( outfile , prop[0] , prop[1] , Op ) ;

  float complex *s = malloc( LCU*NSNS*NCNC * sizeof( float complex ) ) ;
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {

    // by convention the second prop is Anti
    void (*f)( struct spinor *A , const struct spinor B ) ;
    // check if we have gone through a boundary and that the boundaries
    // on the two props are different, if so we flip a sign in prop[1]
    if( ( t >= prop[0].origin[ND-1] ) &&
	( prop[0].bound[ND-1] != prop[1].bound[ND-1] ) ) {
      switch( Op ) {
      case ADD : f = sub_spinors ; break ;
      case SUB : f = add_spinors ; break ;
      case MUL : f = spinmul_atomic_right ; break ;
      }
    } else {
      switch( Op ) {
      case ADD : f = add_spinors ; break ;
      case SUB : f = sub_spinors ; break ;
      case MUL : f = spinmul_atomic_right ; break ;
      }
    }
    
    read_prop( prop[0] , S1 , t ) ;
    read_prop( prop[1] , S2 , t ) ;
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

  if( s != NULL ) {
    free( s ) ;
  }
  
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
