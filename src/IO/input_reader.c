/*
    Copyright 2013 Renwick James Hudspith

    This file (input_reader.c) is part of GLU.

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
   @file input_reader.c
   @brief reader for the input file

   Some of it was taken from my library GLU, some is new
*/

#include "common.h"

#include <errno.h>

// break for the while loops can be made really large
// I just don't like infinite while loops
#define MAX_CONTRACTIONS 20

// tokenize the input file
struct inputs {
  char TOKEN[ GLU_STR_LENGTH ] ;
  char VALUE[ GLU_STR_LENGTH ] ;
} ;

// is this what valgrind dislikes?
static struct inputs *INPUT ;

// counter for the number of tags
static int NTAGS = 0 ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) { return !strcmp( str_1 , str_2 ) ; }

// allocate the input file struct
static void 
pack_inputs( FILE *setup )
{
  INPUT = ( struct inputs* )malloc( MAX_TOKENS * sizeof( struct inputs ) ) ;
  // and put into the buffer
  while( NTAGS++ , fscanf( setup , "%s = %s" , INPUT[ NTAGS ].TOKEN , INPUT[ NTAGS ].VALUE )  != EOF ) { }
  return ;
}

// frees up the struct
static void
unpack_inputs( void )
{
  free( INPUT ) ;
  return ;
}

// prints out the problematic tag
static int
tag_failure( const char *tag )
{
  printf( "[IO] Failure looking for tag %s in input file ... Leaving\n" , tag ) ;
  return FAILURE ;
}

// look for a tag and return the index if found, else return -1
static int
tag_search( const char *tag ) 
{
  int i ;
  for( i = 0 ; i < NTAGS ; i++ ) {
    if( are_equal( INPUT[i].TOKEN , tag ) ) return i ;
  }
  return FAILURE ;
}

// quickly get the configuration number from the input file
static int
confno( void )
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

// get the header type
static header_mode
header_type( void )
{
  const int header_idx = tag_search( "HEADER" ) ;
  if( header_idx == FAILURE ) { return tag_failure( "HEADER" ) ; }

  if( are_equal( INPUT[header_idx].VALUE , "NERSC" ) ) {
    printf( "[IO] Attempting to read a NERSC file \n" ) ;
    return NERSC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "HIREP" ) ) {
    printf( "[IO] Attempting to read a HIREP file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return HIREP_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "MILC" ) ) {
    printf( "[IO] Attempting to read a MILC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return MILC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "SCIDAC" ) ) {
    printf( "[IO] Attempting to read a SCIDAC file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "LIME" ) ) {
    printf( "[IO] Attempting to read an LIME file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    printf( "[IO] WARNING!! NOT CHECKING ANY CHECKSUMS!! \n" ) ;
    return LIME_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_SCIDAC" ) ) {
    printf( "[IO] Attempting to read an ILDG (Scidac) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return ILDG_SCIDAC_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "ILDG_BQCD" ) ) {
    printf( "[IO] Attempting to read an ILDG (BQCD) file \n" ) ;
    printf( "[IO] Using sequence number from input file :: %d \n" ,
	    Latt.flow = confno( ) ) ;
    return ILDG_BQCD_HEADER ;
  } else if( are_equal( INPUT[header_idx].VALUE , "UNIT" ) ) {
    printf( "[IO] Generating a UNIT gauge configuration \n" ) ;
    return UNIT_GAUGE ;
  }
  printf( "[IO] HEADER %s not recognised ... Leaving \n" , 
	  INPUT[header_idx].VALUE ) ;
  return FAILURE ; 
}

// pack the cut_info struct
static int
read_cuts_struct( struct cut_info *CUTINFO )
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

// get the map for our props
static int
get_contraction_map( int *map ,
		     const char *token ,
		     const int nprops ) 
{
  errno = 0 ;
  char *endptr ;
  *map = (int)strtol( token , &endptr , 10 ) ;
  if( token == endptr ) {
    printf( "[IO] contraction mapping expects an integer \n" ) ;
    return FAILURE ;
  }
  if( *map < 0 || *map >= nprops || errno == ERANGE ) {
    printf( "[IO] non-sensical contraction index %d \n" , *map ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// get the various contraction types
static int
get_current_type( current_type *current ,
		  const char *token ,
		  const int meas ) 
{
  if( are_equal( token , "LOCAL_LOCAL" ) ) {
    *current = LOCAL_LOCAL ;
  } else if( are_equal( token , "CONSERVED_LOCAL" ) ) {
    *current = CONSERVED_LOCAL ;
  } else {
    printf( "[IO] I don't understand VPF type %s\n" , token ) ;
    return FAILURE ;
  }
  printf( "[IO] VPF_%d :: We are performing %s contractions \n" , meas , token ) ;
  return SUCCESS ;
}

//
static int
unexpected_NULL( void ) 
{
  printf( "[IO] unexpected NULL in contraction string \n" ) ;
  return FAILURE ;
}

// tokenize each meson
static int
twopoint_tokens( struct meson_info *mesons ,
		 const char *meson_str ,
		 const int nprops ,
		 const int meas_idx ,
		 const char *message ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)meson_str , "," ) ;
  if( get_contraction_map( &( mesons -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( mesons -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] %s_%d :: Contracting prop %d with prop %d \n" , 
	  message , meas_idx , mesons -> map[0] , mesons -> map[1] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( mesons -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%d :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%d :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// tokenize each meson
static int
threepoint_tokens( struct baryon_info *baryons ,
		   const char *baryon_str ,
		   const int nprops ,
		   const int meas_idx ,
		   const char *message ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)baryon_str , "," ) ;
  if( get_contraction_map( &( baryons -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( baryons -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( baryons -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] %s_%d :: Contracting prop %d with prop %d with prop %d \n" , 
	  message , meas_idx , baryons -> map[0] , baryons -> map[1] , baryons -> map[2] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( baryons -> outfile , "%s" , token ) ;
  printf( "[IO] %s_%d :: Contraction file in %s \n" , 
	  message , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] %s_%d :: Unexpected extra contraction info %s \n" ,
	    message , meas_idx , token ) ;
  }
  return SUCCESS ;
}

// baryon contractions of the form prop_idx1,prop_idx2,proptype,outfile
static int
baryon_contractions( struct baryon_info *baryons , 
		     int *nbaryons ,
		     const int nprops ,
		     const GLU_bool first_pass ) 
{
  *nbaryons = 0 ;
  char str[ 32 ] ;
  while( *nbaryons < MAX_CONTRACTIONS ) {
    sprintf( str , "BARYON%d" , *nbaryons ) ;
    const int baryon_idx = tag_search( str ) ;
    if( baryon_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( threepoint_tokens( &baryons[ *nbaryons ] , 
			     INPUT[ baryon_idx ].VALUE , 
			     nprops , *nbaryons , "Baryon" ) 
	  == FAILURE ) return FAILURE ;
    }
    *nbaryons = *nbaryons + 1 ;
  }
  return SUCCESS ;
}


// meson contractions of the form prop_idx1,prop_idx2,proptype,outfile
static int
meson_contractions( struct meson_info *mesons , 
		    int *nmesons ,
		    const int nprops ,
		    const GLU_bool first_pass ) 
{
  *nmesons = 0 ;
  char str[ 32 ] ;
  while( *nmesons < MAX_CONTRACTIONS ) {
    sprintf( str , "MESON%d" , *nmesons ) ;
    const int meson_idx = tag_search( str ) ;
    if( meson_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( twopoint_tokens( &mesons[ *nmesons ] , 
			   INPUT[ meson_idx ].VALUE , 
			   nprops , *nmesons , "Meson" ) 
	  == FAILURE ) return FAILURE ;
    }
    *nmesons = *nmesons + 1 ;
  }
  return SUCCESS ;
}

// tokenize each vpf contraction
static int
VPF_tokens( struct VPF_info *VPF ,
	    const char *VPF_str ,
	    const int nprops ,
	    const int meas_idx ) 
{
  printf( "\n" ) ;

  // starts with the propagator map
  char *token = (char*)strtok( (char*)VPF_str , "," ) ;
  if( get_contraction_map( &( VPF -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( VPF -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] VPF_%d :: Contracting prop %d with prop %d \n" , 
	  meas_idx , VPF -> map[0] , VPF -> map[1] ) ;

  // check for current
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_current_type( &( VPF -> current ) , token , meas_idx ) == FAILURE ) return FAILURE ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( VPF -> outfile , "%s" , token ) ;
  printf( "[IO] VPF_%d :: Contraction file in %s \n" , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] VPF_%d :: Unexpected extra contraction info %s \n" ,
	    meas_idx , token ) ;
  }
  return SUCCESS ;
}

// meson contractions of the form prop_idx1,prop_idx2,proptype,outfile
static int
VPF_contractions( struct VPF_info *VPF , 
		  int *nVPF ,
		  const int nprops ,
		  const GLU_bool first_pass ) 
{
  *nVPF = 0 ;
  char str[ 32 ] ;
  while( *nVPF < MAX_CONTRACTIONS ) {
    sprintf( str , "VPF%d" , *nVPF ) ;
    const int VPF_idx = tag_search( str ) ;
    if( VPF_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( VPF_tokens( &VPF[ *nVPF ] , 
		      INPUT[ VPF_idx ].VALUE , 
		      nprops , *nVPF ) 
	  == FAILURE ) return FAILURE ;
    }
    *nVPF = *nVPF + 1 ;
  }
  return SUCCESS ;
}

// tokenize each meson
static int
matrix_element_tokens( struct WME_info *wme ,
		       const char *WME_str ,
		       const int nprops ,
		       const int meas_idx ) 
{
  printf( "\n" ) ;

  // starts with the contraction indices
  char *token = (char*)strtok( (char*)WME_str , "," ) ;
  if( get_contraction_map( &( wme -> map[0] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[1] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[2] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_contraction_map( &( wme -> map[3] ) , token , nprops ) == FAILURE ) {
    return FAILURE ;
  }
  printf( "[IO] WME_%d :: Contracting prop %d with prop %d with %d with %d \n" , 
	  meas_idx , wme -> map[0] , wme -> map[1] , wme -> map[2] , wme -> map[3] ) ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( wme -> outfile , "%s" , token ) ;
  printf( "[IO] WME_%d :: Contraction file in %s \n" , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] WME_%d :: Unexpected extra contraction info %s \n" ,
	    meas_idx , token ) ;
  }
  return SUCCESS ;
}

// meson contractions of the form prop_idx1,prop_idx2,proptype,outfile
static int
matrix_element_contractions( struct WME_info *wme , 
			     int *nWME ,
			     const int nprops ,
			     const GLU_bool first_pass ) 
{
  *nWME = 0 ;
  char str[ 32 ] ;
  while( *nWME < MAX_CONTRACTIONS ) {
    sprintf( str , "WME%d" , *nWME ) ;
    const int WME_idx = tag_search( str ) ;
    if( WME_idx == FAILURE ) return SUCCESS ;
    if( first_pass == GLU_FALSE ) {
      if( matrix_element_tokens( &wme[ *nWME ] , 
				 INPUT[ WME_idx ].VALUE , 
				 nprops , *nWME ) 
	  == FAILURE ) return FAILURE ;
    }
    *nWME = *nWME + 1 ;
  }
  return SUCCESS ;
}

// get the lattice dimensions from the input_file
static int
get_dims( int *dims )
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
static int
get_props( struct propagator *props ,
	   int *nprops ,
	   const GLU_bool first_pass )
{
  *nprops = 0 ;
  char str[ 32 ] ;
  while( *nprops < 10 ) {
    sprintf( str , "PROP%d" , *nprops ) ;
    const int prop_idx = tag_search( str ) ;
    if( prop_idx == FAILURE ) break ;
    if( first_pass == GLU_FALSE ) {
      char *token ;
      if( ( token = (char*)strtok( INPUT[ prop_idx ].VALUE  , "," ) ) == NULL ) {
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

// little wrapper to free our inputs. I invisage it will grow with time
void
free_inputs( struct input_info inputs ) 
{
  free( inputs.mesons ) ;
  free( inputs.VPF ) ;
  free( inputs.wme ) ;
  return ;
}

// free the propagators, another structure that might get quite big
void
free_props( struct propagator *props , 
	    const int nprops )
{
  int i ;
  for( i = 0 ; i < nprops ; i++ ) {
    fclose( props[i].file ) ;
  }
  free( props ) ;
  return ;
}

// pack the input_info struct with all our contractions
int
get_input_data( struct propagator **prop ,
		struct input_info *inputs ,
		const char *file_name )
{
  // open the input file in here and free it at the bottom
  FILE *infile = fopen( file_name , "r" ) ;
  if( infile == NULL ) {
    printf( "[IO] input file cannot be read ... Leaving\n" ) ;
    return FAILURE ;
  }
  
  //
  int STATUS = SUCCESS ;

  // if we can open the file we push it into a big structure
  pack_inputs( infile ) ;

  // get the header
  Latt.head = header_type( ) ;
  if( Latt.head == FAILURE ) STATUS = FAILURE ;

  // read some props, although strange I wouldn't consider no props being an error
  get_props( *prop , &( inputs -> nprops ) , GLU_TRUE ) ;
  *prop = (struct propagator*)malloc( ( inputs -> nprops ) * sizeof( struct propagator ) ) ;
  if( inputs -> nprops == 0 ) { 
    printf( "[IO] No propagator files specified \n" ) ;
  }
  get_props( *prop , &( inputs -> nprops ) , GLU_FALSE ) ;

  // baryon stuff
  baryon_contractions( inputs -> baryons , &( inputs -> nbaryons ) , inputs -> nprops , GLU_TRUE ) ;
  inputs -> baryons = (struct baryon_info*)malloc( ( inputs -> nbaryons ) * sizeof( struct baryon_info ) ) ;
  if( baryon_contractions( inputs -> baryons , &( inputs -> nbaryons ) , 
			   inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // meson contractions
  meson_contractions( inputs -> mesons , &( inputs -> nmesons ) , inputs -> nprops , GLU_TRUE ) ;
  inputs -> mesons = (struct meson_info*)malloc( ( inputs -> nmesons ) * sizeof( struct meson_info ) ) ;
  if( meson_contractions( inputs -> mesons , &( inputs -> nmesons ) , 
			  inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // vacuum polarisation stuff
  VPF_contractions( inputs -> VPF , &( inputs -> nVPF ) , inputs -> nprops , GLU_TRUE ) ;
  inputs -> VPF = (struct VPF_info*)malloc( ( inputs -> nVPF ) * sizeof( struct VPF_info ) ) ;
  if( VPF_contractions( inputs -> VPF , &( inputs -> nVPF ) , 
			inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // so I think that not having the cuts isn't a problem apart from if we
  // need them, maybe this is a bit dodgy - J
  if( ( inputs -> nVPF ) != 0 && read_cuts_struct( &( inputs -> CUTINFO ) ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // matrix element stuff
  matrix_element_contractions( inputs -> wme , &( inputs -> nWME ) , inputs -> nprops , GLU_TRUE ) ;
  inputs -> wme = (struct WME_info*)malloc( ( inputs -> nWME ) * sizeof( struct WME_info ) ) ;
  if( matrix_element_contractions( inputs -> wme , &( inputs -> nWME ) , 
				   inputs -> nprops , GLU_FALSE ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // read the dimensions from the input file
  if( get_dims( inputs -> dims ) == FAILURE ) {
    STATUS = FAILURE ;
  }

  // close the file and deallocate the buffer
  fclose( infile ) ;
  unpack_inputs( ) ;

  // if we have hit ANY problem we return GLU_FAILURE this causes it to exit
  return STATUS ;
}

