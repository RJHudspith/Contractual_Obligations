/**
   @file corr.c
   @brief mainfile
 */
#include "common.h"          // one header to rule them all

#include "conserved_local.h" // conserved local currents
#include "geometry.h"        // init_geom and init_navig
#include "GLU_timer.h"       // sys/time.h wrapper
#include "io.h"              // file IO stuff
#include "input_reader.h"    // input file reader
#include "mesons.h"          // meson contractions
#include "plaqs_links.h"     // plaquettes and links of gauge field
#include "readers.h"         // gauge config reader
#include "read_config.h"     // read a gauge configuration file
#include "read_headers.h"    // read the header file in gauge config
#include "read_propheader.h" // read the propagator file header
#include "brutal_mesons.h"   // brutal mesons

// lattice information holds dimensions and other stuff
// to be taken from the gauge configuration file OR the input file
struct latt_info Latt ;

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

  // we should loop this section to read in many files
  int nprops = 4 , dims[ ND ] ;
  char prop_files[ nprops ][ GLU_STR_LENGTH ] ;
  if( get_input_data( prop_files , &nprops , dims , argv[INFILE] ) == FAILURE ) {
    return FAILURE ;
  }

  // geometry has to come from the input file
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = dims[ mu ] ;
  }
  init_geom( ) ;

  // at the moment we read in the whole gauge field
  // maybe that will change I am not sure, we could probably live with it
  // like this, the gauge field is a big beast but not as big as a prop
  // and we can probably free it if needs be
  struct site *lat = NULL ;
  struct head_data HEAD_DATA ;
  if( MODE == GAUGE_AND_PROPS ) {
    lat = read_gauge_file( &HEAD_DATA , argv[GAUGE_FILE] ) ;
    // check input file versus gauge field dimensions? Probably should
  } 

  // open up some propagator files
  FILE *fprops[ nprops ] ;
  int i = 0 ;
  for( i = 0 ; i < nprops ; i++ ) {
    // open and check all the files
    if( ( fprops[i] = fopen( prop_files[i] , "r" ) ) == NULL ) {
      printf( "[IO] Propagator file %s empty! Leaving \n" , prop_files[i] ) ;
      if( MODE == GAUGE_AND_PROPS ) free( lat ) ;
      return FAILURE ;
    }
    //
    if( read_check_header( fprops[i] ) == FAILURE ) {
      if( MODE == GAUGE_AND_PROPS ) free( lat ) ;
      return FAILURE ;
    }
  }

  start_timer( ) ;

  // want to switch on these or call a wrapper
  //single_mesons( fprops[0] , CHIRAL ) ;
  //single_mesons_bruteforce( fprops[0] , CHIRAL ) ;
  //hheavy_mesons( fprops[0] , NREL ) ;
  //double_mesons( fprops[0] , CHIRAL_TO_NREL , fprops[1] , NREL ) ;
  //double_mesons( fprops[2] , NREL , fprops[3] , NREL ) ;

  if( lat != NULL ) {
    conserved_local( fprops[0] , CHIRAL , lat ) ;
  }

  //wall_mesons( fprops[0] , CHIRAL ) ;

  print_time( ) ;

  // we will have to move this around only place where this is freed
  if( MODE == GAUGE_AND_PROPS ) {
    free( lat ) ;
  }

  for( i = 0 ; i < nprops ; i++ ) {
    fclose( fprops[i] ) ;
  }

  return SUCCESS ;
}
