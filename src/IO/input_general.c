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
size_t
confno( const struct inputs *INPUT )
{
  errno = 0 ;
  char *endptr ;
  const int conf_idx = tag_search( "CONFNO" ) ;
  const size_t num = (size_t)strtol( INPUT[conf_idx].VALUE , 
				     &endptr , 10 ) ; 
  if( endptr == INPUT[conf_idx].VALUE || errno == ERANGE ) {
    return 0 ;
  }
  // should also check that it is a sensible int as we do a cast ...
  return num ;
}

// get the DFT information
static int
get_DFT( double *proto_mom ,
	 double *thetas ,
	 size_t *nalphas ,
	 const struct inputs *INPUT )
{
  errno = 0 ;
  char *token , *endptr , str[ 32 ] ;
  const int dft_idx = tag_search( "DFT" ) ;
  if( dft_idx == FAILURE ) { 
    return tag_failure( str ) ; 
  }
  token = strtok( (char*)INPUT[dft_idx].VALUE , "," ) ;
  size_t mu = 0 ;
  proto_mom[ mu ] = strtod( token , &endptr ) ;
  mu++ ;
  while( ( token = strtok( NULL , "," ) ) != NULL ) {
    if( mu < ND-1 ) {
      proto_mom[ mu ] = strtod( token , &endptr ) ;
    } else {
      thetas[ mu - (ND-1) ] = strtod( token , &endptr ) ;
    }
    mu++ ;
  }
  *nalphas = mu - (ND-1);
  printf( "[IO] Nalphas %zu \n" , *nalphas ) ;
  return SUCCESS ;
}

// get the lattice dimensions from the input_file
int
get_dims( size_t *dims , 
	  const struct inputs *INPUT )
{
  errno = 0 ;
  char *endptr , str[ 32 ] ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    sprintf( str , "DIMS%zu" , mu ) ;
    const int dims_idx = tag_search( str ) ;
    if( dims_idx == FAILURE ) { 
      return tag_failure( str ) ; 
    }
    dims[ mu ] = (size_t)strtol( INPUT[ dims_idx ].VALUE , &endptr , 10 ) ;
    if( INPUT[ dims_idx ].VALUE == endptr ||
	errno == ERANGE ) {
      fprintf( stderr , "[IO] non-sensical dimension value %zu \n" , 
	       dims[ mu ] ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// get prop tags
int
get_props( struct propagator *props ,
	   size_t *nprops ,
	   const struct inputs *INPUT ,
	   const GLU_bool first_pass )
{
  *nprops = 0 ;
  char str[ 32 ] ;
  while( *nprops < MAX_CONTRACTIONS ) {
    sprintf( str , "PROP%zu" , *nprops ) ;
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
	fprintf( stderr , "[IO] propfile %s not found \n" , token ) ;
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
    fprintf( stdout , "[IO] Attempting to read a NERSC file \n" ) ;
    return NERSC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "CERN" ) ) {
    fprintf( stdout , "[IO] Attempting to read a CERN file \n" ) ;
    return CERN_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "HIREP" ) ) {
    fprintf( stdout , "[IO] Attempting to read a HIREP file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    return HIREP_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "MILC" ) ) {
    fprintf( stdout , "[IO] Attempting to read a MILC file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    return MILC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "SCIDAC" ) ) {
    fprintf( stdout , "[IO] Attempting to read a SCIDAC file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    return SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "LIME" ) ) {
    fprintf( stdout , "[IO] Attempting to read an LIME file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    printf( "[IO] WARNING!! NOT CHECKING ANY CHECKSUMS!! \n" ) ;
    return LIME_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_SCIDAC" ) ) {
    fprintf( stdout , "[IO] Attempting to read an ILDG (Scidac) file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    return ILDG_SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_BQCD" ) ) {
    fprintf( stdout , "[IO] Attempting to read an ILDG (BQCD) file \n" ) ;
    fprintf( stdout , "[IO] Using sequence number from input file :: %zu \n" ,
	     Latt.flow = confno( INPUT ) ) ;
    return ILDG_BQCD_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "UNIT" ) ) {
    fprintf( stdout , "[IO] Generating a UNIT gauge configuration \n" ) ;
    return UNIT_GAUGE ;
  }
  fprintf( stderr , "[IO] HEADER %s not recognised ... Leaving \n" , 
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
  } else if ( are_equal( INPUT[momcut_idx].VALUE , "PSQ_CUT" ) ) {
    CUTINFO -> type = PSQ_CUT ; 
  } else if ( are_equal( INPUT[momcut_idx].VALUE , "CYLINDER_CUT" ) ) {
    CUTINFO -> type = CYLINDER_CUT ; 
  } else {
    fprintf( stderr , "[IO] Unrecognised type [%s] \n" , 
	     INPUT[momcut_idx].VALUE ) ; 
    fprintf( stderr , "[IO] Defaulting to PSQ_CUT \n" ) ; 
    CUTINFO -> type = PSQ_CUT ; 
  }
  // minmom, maxmom angle and cylinder width
  const int maxmom_idx = tag_search( "MAXMOM" ) ;
  if( maxmom_idx == FAILURE ) { return tag_failure( "MAXMOM" ) ; }
  errno = 0 ;
  CUTINFO -> max_mom = (int)strtol( INPUT[maxmom_idx].VALUE , &endptr , 10 ) ;
  if( endptr == INPUT[maxmom_idx].VALUE || errno == ERANGE ) {
    printf( "[IO] non-sensical maximum momentum %zu \n" , CUTINFO -> max_mom ) ;
    return FAILURE ;
  }  
  const int cyl_idx = tag_search( "CYL_WIDTH" ) ;
  if( cyl_idx == FAILURE ) { return tag_failure( "CYL_WIDTH" ) ; }
  errno = 0 ;
  CUTINFO -> cyl_width = (double)strtod( INPUT[cyl_idx].VALUE , &endptr ) ;
  if( endptr == INPUT[cyl_idx].VALUE || errno == ERANGE ||
      CUTINFO -> cyl_width < 0.0 ) {
    fprintf( stderr , "[IO] non-sensical cylinder width %f \n" , 
	     CUTINFO -> cyl_width ) ;
    return FAILURE ;
  }
  // DFT information
  if( get_DFT( CUTINFO -> proto_mom , CUTINFO -> thetas ,
	       &CUTINFO -> Nalphas , INPUT ) == FAILURE ) {
    fprintf( stdout , "[IO] Not performing DFT\n" ) ;
    CUTINFO -> Nalphas = 0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      CUTINFO -> proto_mom[mu] = 0.0 ;
    }
  }
  // config space
  const int cspace_idx = tag_search( "CONFIGSPACE" ) ;
  if( are_equal( INPUT[cspace_idx].VALUE , "TRUE" ) ) {
    fprintf( stdout , "[IO] computing configuration space correlators\n" ) ;
    CUTINFO -> configspace = GLU_TRUE ;
  } else {
    fprintf( stdout , "[IO] computing momentum space correlators\n" ) ;
    CUTINFO -> configspace = GLU_FALSE ;
  }
  // maximum R2
  const int maxr2_idx = tag_search( "MAXR2" ) ;
  errno = 0 ;
  CUTINFO -> max_r2 = (int)strtol( INPUT[maxr2_idx].VALUE , &endptr , 10 ) ;
  if( endptr == INPUT[maxr2_idx].VALUE || errno == ERANGE ) {
    printf( "[IO] non-sensical maximum r^2 %zu \n" , CUTINFO -> max_r2 ) ;
    return FAILURE ;
  }

  return SUCCESS ;
}
