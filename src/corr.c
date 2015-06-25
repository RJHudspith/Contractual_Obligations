/**
   @file corr.c
   @brief mainfile
 */
#include "common.h"          // one header to rule them all

#include "geometry.h"        // init_geom and init_navig
#include "GLU_timer.h"       // sys/time.h wrapper
#include "input_reader.h"    // input file readers
#include "read_config.h"     // read a gauge configuration file
#include "read_propheader.h" // read the propagator file header

#include "polyakov.h"

#include "wrap_baryons.h"    // Baryon contraction wrapper
#include "wrap_mesons.h"     // Meson contraction wrapper
#include "wrap_tetras.h"     // Tetraquark contraction wrapper
#include "wrap_VPF.h"        // VPF contraction wrapper
#include "wrap_WME.h"        // VPF contraction wrapper

// lattice information holds dimensions and other stuff
// to be taken from the gauge configuration file OR the input file
struct latt_info Latt ;

// gauge field is global for now
struct site *lat = NULL ;

// enumeration for the arguments to our binary
enum{ INFILE = 2 , GAUGE_FILE = 4 } ;

// enumerated mode type
typedef enum{ PROPAGATORS_ONLY , GAUGE_AND_PROPS } corrmode ;

// our main function
int 
main( const int argc, 
      const char *argv[] ) 
{
  // defines what mode we are in
  corrmode MODE = GAUGE_AND_PROPS ;
  if( argc == ( INFILE + 1 ) ) {
    MODE = PROPAGATORS_ONLY ;
  } else if( argc != ( GAUGE_FILE + 1 ) ) {
    printf( "USAGE :: ./CORR -i {input_file} -g {gauge_field} \n" ) ;
    return FAILURE ;
  }

  // read the inputs in the input file
  struct propagator *prop ; // gets allocated in input data
  struct input_info inputs ;
  if( get_input_data( &prop , &inputs ,
		      argv[INFILE] ) == FAILURE ) {
    free_props( prop , inputs.nprops ) ;
    free_inputs( inputs ) ; 
    return FAILURE ;
  }

  // geometry has to come from the input file
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = inputs.dims[ mu ] ;
  }
  init_geom( ) ;

  // my gauge field code requires us to read in the whole config
  struct head_data HEAD_DATA ;
  if( MODE == GAUGE_AND_PROPS ) {
    if( ( lat = read_gauge_file( &HEAD_DATA , argv[GAUGE_FILE] ,
				 inputs.dims ) ) == NULL ) {
      free_props( prop , inputs.nprops ) ;
      free_inputs( inputs ) ; 
      return FAILURE ;
    }
  }

  // open up some propagator files and parse the
  // header checking the geometry
  if( read_propheaders( prop , inputs ) == FAILURE ) {
    if( MODE == GAUGE_AND_PROPS ) free( lat ) ;
    free_props( prop , inputs.nprops ) ;
    free_inputs( inputs ) ;
    return FAILURE ;
  }

  start_timer( ) ;

  // baryon contraction code
  if( contract_baryons( prop , inputs.baryons , inputs.CUTINFO ,
			inputs.nbaryons ) 
      == FAILURE ) {
    goto FREES ; // do not pass GO, do not collect £200
  }

  // calls a wrapper that has some logic tests
  if( contract_mesons( prop , inputs.mesons , inputs.CUTINFO ,
		       inputs.nmesons ) == FAILURE ) {
    goto FREES ; // do not pass GO, do not collect £200
  }

  // tetraquark contraction code
  if( contract_tetras( prop , inputs.tetras , inputs.CUTINFO , 
		       inputs.ntetras ) == FAILURE ) {
    goto FREES ; // do not pass GO, do not collect £200
  }

  // if we don't have a gauge field we can't do conserved-local
  if( lat != NULL ) {
    if( contract_VPF( prop , lat , inputs.VPF ,
		      inputs.nVPF , inputs.CUTINFO ) == FAILURE ) {
      goto FREES ; // do not pass GO, do not collect £200
    }
  } 

  // WME contraction, props have to be wall source
  if( contract_WME( prop , inputs.wme , 
		    inputs.nWME ) == FAILURE ) {
    goto FREES ; // do not pass GO, do not collect £200
  }

 FREES :
  // we will have to move this around only place where this is freed
  free( lat ) ;

  // free the propagators
  free_props( prop , inputs.nprops ) ;

  // free our contraction tables
  free_inputs( inputs ) ; 

  return SUCCESS ;
}
