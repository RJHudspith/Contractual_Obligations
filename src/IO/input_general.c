/**
   @file input_general.c
   @brief general input file stuff like lattice geometry and config headers
 */
#include "common.h"

#include <errno.h>
#include "input_reader.h"

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// quickly get the configuration number from the input file
int
confno( const struct inputs *INPUT )
{
  errno = 0 ;
  char *endptr ;
  const int conf_idx = tag_search( "CONFNO" ) ;
  const int num = (int)strtol( INPUT[conf_idx].VALUE , 
			       &endptr , 10 ) ; 
  if( endptr == INPUT[conf_idx].VALUE || errno == ERANGE ) {
    return FAILURE ;
  }
  // should also check that it is a sensible int as we do a cast ...
  return num ;
}

// get the lattice dimensions from the input_file
int
get_dims( int *dims , 
	  const struct inputs *INPUT )
{
  errno = 0 ;
  char *endptr , str[ 32 ] ;
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    sprintf( str , "DIMS%d" , mu ) ;
    const int dims_idx = tag_search( str ) ;
    if( dims_idx == FAILURE ) { 
      return tag_failure( str ) ; 
    }
    dims[ mu ] = (int)strtol( INPUT[ dims_idx ].VALUE , &endptr , 10 ) ;
    if( dims[ mu ] < 0 || INPUT[ dims_idx ].VALUE == endptr ||
	errno == ERANGE ) {
      printf( "[IO] non-sensical dimension value %d \n" , dims[ mu ] ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// get prop tags
int
get_props( struct propagator *props ,
	   int *nprops ,
	   const struct inputs *INPUT ,
	   const GLU_bool first_pass )
{
  *nprops = 0 ;
  char str[ 32 ] ;
  while( *nprops < MAX_CONTRACTIONS ) {
    sprintf( str , "PROP%d" , *nprops ) ;
    const int prop_idx = tag_search( str ) ;
    if( prop_idx == FAILURE ) break ;
    if( first_pass == GLU_FALSE ) {
      char *token ;
      if( ( token = (char*)strtok( (char*)INPUT[ prop_idx ].VALUE  , "," ) )== NULL ) {
	return unexpected_NULL( ) ;
      }
      // open file
      props[ *nprops ].file = fopen( token , "rb" ) ;
      if( props[ *nprops ].file == NULL ) {
	printf( "[IO] propfile %s not found \n" , token ) ;
	return FAILURE ;
      }
    }
    *nprops = *nprops + 1 ;
  }
  return SUCCESS ;
}

// get the header type
header_mode
header_type( const struct inputs *INPUT )
{
  const int header_idx = tag_search( "HEADER" ) ;
  if( header_idx == FAILURE ) { 
    tag_failure( "HEADER" ) ; 
    return UNSUPPORTED ;
  }
  if( are_equal( INPUT[header_idx].VALUE , "NERSC" ) ) {
    printf( "[IO] Attempting to read a NERSC file \n" ) ;
    return NERSC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "HIREP" ) ) {
    printf( "[IO] Attempting to read a HIREP file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    return HIREP_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "MILC" ) ) {
    printf( "[IO] Attempting to read a MILC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    return MILC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "SCIDAC" ) ) {
    printf( "[IO] Attempting to read a SCIDAC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    return SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "LIME" ) ) {
    printf( "[IO] Attempting to read an LIME file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    printf( "[IO] WARNING!! NOT CHECKING ANY CHECKSUMS!! \n" ) ;
    return LIME_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_SCIDAC" ) ) {
    printf( "[IO] Attempting to read an ILDG (Scidac) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    return ILDG_SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_BQCD" ) ) {
    printf( "[IO] Attempting to read an ILDG (BQCD) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( INPUT ) ) ;
    return ILDG_BQCD_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "UNIT" ) ) {
    printf( "[IO] Generating a UNIT gauge configuration \n" ) ;
    return UNIT_GAUGE ;
  }
  printf( "[IO] HEADER %s not recognised ... Leaving \n" , 
	  INPUT[header_idx].VALUE ) ;
  return UNSUPPORTED ; 
}

// pack the cut_info struct
int
read_cuts_struct( struct cut_info *CUTINFO ,
		  const struct inputs *INPUT )
{
  char *endptr ; // end pointer collection for strto*
  // momentum space cut def
  const int momcut_idx = tag_search( "MOM_CUT" ) ;
  if( momcut_idx == FAILURE ) { return tag_failure( "MOM_CUT" ) ; }
  if ( are_equal( INPUT[momcut_idx].VALUE , "HYPERCUBIC_CUT" ) ) {
    CUTINFO -> type = HYPERCUBIC_CUT ; 
  } else if ( are_equal( INPUT[momcut_idx].VALUE , "SPHERICAL_CUT" ) ) {
    CUTINFO -> type = PSQ_CUT ; 
  } else if ( are_equal( INPUT[momcut_idx].VALUE , "CYLINDER_CUT" ) ) {
    CUTINFO -> type = CYLINDER_CUT ; 
  } else {
    printf( "[IO] Unrecognised type [%s] \n" , INPUT[momcut_idx].VALUE ) ; 
    printf( "[IO] Defaulting to SPHERICAL_CUT \n" ) ; 
    CUTINFO -> type = PSQ_CUT ; 
  }
  // minmom, maxmom angle and cylinder width
  const int maxmom_idx = tag_search( "MAXMOM" ) ;
  if( maxmom_idx == FAILURE ) { return tag_failure( "MAXMOM" ) ; }
  errno = 0 ;
  CUTINFO -> max_mom = (int)strtol( INPUT[maxmom_idx].VALUE , &endptr , 10 ) ;
  if( endptr == INPUT[maxmom_idx].VALUE || errno == ERANGE || 
      CUTINFO -> max_mom < 1 ) {
    printf( "[IO] non-sensical maximum momentum %d \n" , CUTINFO -> max_mom ) ;
    return FAILURE ;
  }  
  const int cyl_idx = tag_search( "CYL_WIDTH" ) ;
  if( cyl_idx == FAILURE ) { return tag_failure( "CYL_WIDTH" ) ; }
  errno = 0 ;
  CUTINFO -> cyl_width = (double)strtod( INPUT[cyl_idx].VALUE , &endptr ) ;
  if( endptr == INPUT[cyl_idx].VALUE || errno == ERANGE ||
      CUTINFO -> cyl_width < 0.0 ) {
    printf( "[IO] non-sensical cylinder width %f \n" , CUTINFO -> cyl_width ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}
