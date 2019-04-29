/*
    Copyright 2013 Renwick James Hudspith

    This file (read_config.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file read_config.c 
   @brief gets the information about our configs from the header
 */
#include "common.h"

#include "CERN.h"          // read_CLS_field()
#include "corr_malloc.h"   // align to the correct boundary
#include "geometry.h"      // init_navig is here
#include "GLU_timer.h"     // print_time()
#include "HIREP.h"         // read a hirep gauge config
#include "plaqs_links.h"   // average plaquette, link traces
#include "readers.h"       // read config files
#include "read_headers.h"  // header readers
#include "Scidac.h"        // Scidac header reading

// identity matrix
static void
identity( double complex a[ NCNC ] )
{
#if NC == 3
  a[0] = 1.0 ; a[1] = 0.0 ; a[2] = 0.0 ;
  a[3] = 0.0 ; a[4] = 1.0 ; a[5] = 0.0 ;
  a[6] = 0.0 ; a[7] = 0.0 ; a[8] = 1.0 ;
#else
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    a[i] = 0.0 ;
  }
  for( i = 0 ; i < NC ; i++ ) {
    a[ i*(NC+1) ] = 1.0 ;
  }
#endif
  return ;
}

// All ones 
static void
unit_gauge( struct site *__restrict lat )
{
  size_t i ; 
  fprintf( stdout , "\n[UNIT] Creating identity SU(%d) lattice fields \n" , 
	   NC ) ; 
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) { 
      identity( lat[i].O[mu] ) ; 
    } 
  }    
  return ;
}

// comparison between checksums
static int 
check_sums( const double plaq , 
	    const double tr ,
	    const uint32_t chksum ,
	    const struct head_data HEAD_DATA )
{
  // these headers I use the same check, the %29 one
  if( Latt.head == MILC_HEADER || Latt.head == ILDG_SCIDAC_HEADER || 
      Latt.head == SCIDAC_HEADER ) {
    // do a check only on the sum29
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums [Calc] %x [Read] %x \n\n" , 
	      chksum , HEAD_DATA.checksum ) ; 
      return FAILURE ;
    }
  } else if( Latt.head == ILDG_BQCD_HEADER ) {
    // BQCD provides a CRC and a plaquette
    double TTOL = 0.0 ;
    if( HEAD_DATA.precision == FLOAT_PREC ) {
      TTOL = 1E-6 ;
    } else {
      TTOL = 1E-14 ;
    }
    if( fabs( HEAD_DATA.plaquette - plaq ) > TTOL ) {
      fprintf( stderr , "\n[IO] Unequal Plaquettes %e %e \n\n" , 
	       plaq , HEAD_DATA.plaquette ) ; 
      return FAILURE ;
    }
    // and the CRC of the binary data
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	       chksum , HEAD_DATA.checksum ) ; 
      return FAILURE ;
    }
  } else if( Latt.head == HIREP_HEADER ) {
    // only check available is the plaquette
    if( fabs( plaq - HEAD_DATA.plaquette ) > PREC_TOL ) {
      fprintf( stderr , "[IO] HIREP header Plaquette Mismatch"
	       "%e vs %e ... Leaving \n" ,
	       plaq , HEAD_DATA.plaquette ) ;
      return FAILURE ;
    }
  } else if( Latt.head == NERSC_HEADER ) {
    enum{ DISASTER = 3 } ;
    // it is disastrous if all three checks fail ...
    int error = 0 ; 
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	       chksum , HEAD_DATA.checksum ) ; 
      error ++ ; 
      // TOL is defined as 10^-6
    } if( fabs( plaq - HEAD_DATA.plaquette ) > PLAQ_AND_TRACE_TOL ) {
      fprintf( stderr , "\n[IO] Unequal Plaquettes Calc %f || Read %f \n\n" , 
	       plaq , HEAD_DATA.plaquette ) ; 
      error ++ ; 
    } if( fabs( tr - HEAD_DATA.trace) > PLAQ_AND_TRACE_TOL ) {
      fprintf( stderr , "\n[IO] Unequal Link_Traces Calc %1.8f "
	       "|| Read %1.8f \n\n" , 
	       tr , HEAD_DATA.trace ) ; 
      error ++ ; 
    }
    // pretty printing
    fprintf( stderr , "[IO] Header     Trace :: %f           || Plaq :: %f \n" ,
	     HEAD_DATA.trace , HEAD_DATA.plaquette ) ; 
    // if everything is wrong we leave
    if( error == DISASTER ) {
      fprintf( stderr , "[IO] NONE of the NERSC Checksums match,"
	       " this is a problem .. Leaving \n") ; 
      return FAILURE ;
    }
  } else if( Latt.head == CERN_HEADER ) {
    fprintf( stderr , "[IO] Header Plaq :: %1.15f \n" , HEAD_DATA.plaquette ) ; 
    if( fabs( plaq - HEAD_DATA.plaquette ) > PLAQ_AND_TRACE_TOL ) {
      fprintf( stderr , "\n[IO] Unequal Plaquettes Calc %f || Read %f \n\n" , 
	       plaq , HEAD_DATA.plaquette ) ;
      return FAILURE ;
    }
  }
  fprintf( stdout , "[IO] Calculated Trace :: %1.15f  || Plaq :: %1.15f \n" , 
	   tr , plaq ) ; 
  return SUCCESS ; // may only be partially successful but I am optimistic
}

// defined in the autoconf procedure
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
  #define DANGEROUS // if we are less than the page size we plough on

  #ifdef GLU_BGQ // BGQ says it has twice as much RAM as it (usually for edinburgh) has ...
  static const double GB = 0.5 * 9.3132257461547852e-10 ; // the mysterious factor is 1.0/(1024*1024*1024)
  #else
  static const double GB = 9.313225746147852e-10 ;
  #endif
static double page_size , MemTotal , MemFree , gauge_fields , gtrans_mats , 
              lat_element , sublat_element ;
static int IO = -1 ;

typedef enum
  { FAST ,
    MODERATE ,
    SLOW } GLU_speed ;

// call to the OS to request what memory size we have
static void
check_mem( void ) 
{
  const double GLU_size = (double)sizeof( double complex ) ;
  page_size = (double)sysconf( _SC_PAGESIZE ) ; 
  MemTotal = sysconf(_SC_PHYS_PAGES) * page_size ; 
  MemFree = (double)sysconf( _SC_AVPHYS_PAGES ) * page_size ; 
  gauge_fields = GLU_size * (double)LVOLUME * ( NCNC * ND ) * GB + 
    (double)sizeof( int ) * (double)LVOLUME * ( 2 * ND ) * GB ; // integer neighbours
  gtrans_mats = (double)LVOLUME * GLU_size * ( NCNC ) * GB ;
  // lattice sites in GB
  lat_element = GLU_size * (double)LVOLUME * GB ;
  // subcube in GB
  sublat_element = GLU_size * (double)LCU * GB ;
  return ;
}
// and print out the standard message
static void
print_info( const char *info ,
	    const double free_memory ,
	    const double expensive_more ,
	    const double expensive_less )
{
  fprintf( stdout , "[%s] Page size is :: %f, Total malloc "
	   "[More] %f  [Less] %f \n" , 
	   info , ( double )free_memory , expensive_more , expensive_less ) ; 
  return ;
}
#endif

static short int 
have_memory_readers_writers( const GLU_output config_type )
{
  int CONFIG = FAST ;
#ifdef HAVE_UNISTD_H
  if( IO == -1 ) {
    check_mem( ) ;
    
    double expensive_more , expensive_less ;
    
    #ifdef DANGEROUS
    const double free_memory = MemTotal * GB ;
    #else
    const double free_memory = MemFree * GB ;
    #endif

    // check to make sure we can fit gauge field and the malloc
    IO = FAST ;
    switch( config_type ) {
    case OUTPUT_SMALL :
      expensive_more = gauge_fields + (double)sizeof( double ) * ( (double)LVOLUME * ( ND * ( NCNC - 1 ) ) ) * GB ;
      expensive_less = gauge_fields + (double)sizeof( double ) * ( (double)( ND * ( NCNC - 1 ) ) ) * GB ;
      print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
      if( expensive_more < free_memory ) {
	  return FAST ;
      } else { 
	IO = SLOW ;
	return SLOW ;
      }
    case OUTPUT_GAUGE :
      expensive_more = gauge_fields + (double)sizeof( double complex ) * ( (double)LVOLUME * ( ND * NC * ( NC - 1 ) ) ) * GB ;
      expensive_less = gauge_fields + (double)sizeof( double complex ) * ( (double)( ND * NC * ( NC - 1 ) ) ) * GB ;
      print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
      if( expensive_more < free_memory ) {
	return FAST ;
      } else { 
	IO = SLOW ;
	return SLOW ;
      }
    case OUTPUT_MILC :
    case OUTPUT_SCIDAC :
    case OUTPUT_ILDG :
    case OUTPUT_NCxNC :
      expensive_more = 2.0 * gauge_fields ; 
      expensive_less = gauge_fields + (double) sizeof( double complex ) * ( (double)( ND * NCNC ) ) * GB ;
      print_info( "IO" , free_memory , expensive_more , expensive_less ) ;
      if( expensive_more < free_memory ) {
	return FAST ;
      } else { 
	IO = SLOW ;
	return SLOW ;
      }
    case OUTPUT_HIREP :
      return FAST ;
    default :
      return SLOW ;
    } 
  } else { CONFIG = IO ; }
#endif
  return CONFIG ;
} 

// things that check the checksums //
short int
checks( struct site *__restrict lat , 
	uint32_t chksum ,
	struct head_data HEAD_DATA )
{
  const double plaq = av_plaquette( lat ) ; 
  const double tr = links( lat ) ; 
  return check_sums( plaq , tr , chksum , HEAD_DATA ) ;
}

// we wrap this one to our new reader ... 
uint32_t
get_config_SUNC( FILE *__restrict CONFIG , 
		 struct site *__restrict lat ,
		 const struct head_data HEAD_DATA )
{
  uint32_t chksum = 0.0 ;
  switch( Latt.head ) {
  case CERN_HEADER :
    read_CLS_field( lat , CONFIG , &chksum ) ;
    return chksum ;
  case LIME_HEADER : // is the same but doesn't care about the checksums
  case ILDG_BQCD_HEADER : // basically all the same NERSC NCxNC
  case ILDG_SCIDAC_HEADER : // ILDG and SCIDAC
  case SCIDAC_HEADER : // Scidac's binary data is compliant
  case MILC_HEADER : // MILC's binary data is the same
  case NERSC_HEADER :
    if( HEAD_DATA.config_type == OUTPUT_SMALL ||
	HEAD_DATA.config_type == OUTPUT_GAUGE ||
	HEAD_DATA.config_type == OUTPUT_NCxNC ) {
      if( have_memory_readers_writers( HEAD_DATA.config_type ) == FAST ) {
	fprintf( stdout , "[IO] using the faster, more wasteful gauge"
		 " configuration reader\n" ) ;
	return lattice_reader_suNC( lat , CONFIG , HEAD_DATA ) ;
      } else {
	fprintf( stdout , "[IO] using the slow site-by-site gauge"
		 " configuration reader\n" ) ;
	return lattice_reader_suNC_cheaper( lat , CONFIG , HEAD_DATA ) ;
      }
    } break ;
  case HIREP_HEADER :
    return read_gauge_field( lat , CONFIG ) ;
  case UNIT_GAUGE :
    unit_gauge( lat ) ;
    return SUCCESS ;
  default : 
    fprintf( stderr , "[IO] Unrecognised HEADER type .. Leaving \n" ) ;
    return FAILURE ;
  }
  return FAILURE ;
}

// read a file, has to be out of order because it is called by the others
struct site*
read_gauge_file( struct head_data *HEAD_DATA , 
		 const char *config_in ,
		 const size_t refdims[ ND ] )
{
  /// here we include the usual stuff look at header for global stuff
  // open our configuration
  FILE *infile = fopen( config_in , "r" ) ;
  if ( infile == NULL ) {
    fprintf( stderr , "[IO] error opening file :: %s\n" , config_in ) ;
    return NULL ;
  }
 
  // initialise the configuration number to zero
  struct head_data tmp ;
  if( read_header( infile , &tmp , GLU_TRUE ) == FAILURE ) {
    fprintf( stderr , "[IO] Header reading failure \n" ) ;
    fclose( infile ) ;
    return NULL ;
  } 

  // initialise geometry so that we can use LVOLUME and stuff
  init_geom( ) ;

  // read_gauge_file overwrites Latt.dims, check against input file
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( Latt.dims[ mu ] != refdims[ mu ] ) {
      fprintf( stderr , "[IO] gauge_field and input file"
	       " dimensions mismatch %zu != %zu \n" 
	       , Latt.dims[ mu ] , refdims[ mu ] ) ;
      fclose( infile ) ;
      return NULL ;
    }
  }

  start_timer() ;

  // malloc our gauge field and initialise our lattice geometry
  struct site *lat = NULL ;
  if( corr_malloc( (void**)&lat , ALIGNMENT ,
		   LVOLUME * sizeof( struct site ) ) != 0 ) {
    fprintf( stderr , "[IO] Gauge field allocation failure ... Leaving \n" ) ;
    fclose( infile ) ;
    free( lat ) ;
    return NULL ;
  }
  init_navig( lat ) ;

  const int check = get_config_SUNC( infile , lat , tmp ) ;
  // read in the configuration ...  
  if( check == FAILURE ) {
    fprintf( stderr , "[IO] File read error ... Leaving \n" ) ;
    fclose( infile ) ;
    free( lat ) ;
    return NULL ;
  }

  // look at scidac header again to get the checksums
  // this is taken from the bottom of the file
  if( Latt.head == SCIDAC_HEADER || Latt.head == ILDG_SCIDAC_HEADER ||
      Latt.head == ILDG_BQCD_HEADER ) {
    get_header_data_SCIDAC( infile , &tmp ) ;
  }

  // have a look at some available checks
  if( checks( lat , check , tmp ) == FAILURE ) { 
    fclose( infile ) ;
    free( lat ) ;
    return NULL ; 
  }

  print_time() ;

  // set the header info
  *HEAD_DATA = tmp ;

  fclose( infile ) ;

  // and finally set the header data into a constant struct
  return lat ;
}
