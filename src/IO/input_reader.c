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
  return atoi( INPUT[tag_search( "CONFNO" )].VALUE ) ;
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
  CUTINFO -> max_mom = atoi( INPUT[maxmom_idx].VALUE ) ;
  const int cyl_idx = tag_search( "CYL_WIDTH" ) ;
  if( cyl_idx == FAILURE ) { return tag_failure( "CYL_WIDTH" ) ; }
  CUTINFO -> cyl_width = atof( INPUT[cyl_idx].VALUE ) ;
  return SUCCESS ;
}

// get the map for our props
static int
get_contraction_map( int *map ,
		     const char *token ,
		     const int nprops ) 
{
  *map = (int)atoi( token ) ;
  if( *map < 0 || *map >= nprops ) {
    printf( "[IO] non-sensical contraction index %d \n" , *map ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// get the various contraction types
static int
get_current_type( current_type *current ,
		  const char *token ) 
{
  if( are_equal( token , "LOCAL_LOCAL" ) ) {
    *current = LOCAL_LOCAL ;
  } else if( are_equal( token , "CONSERVED_LOCAL" ) ) {
    *current = CONSERVED_LOCAL ;
  } else {
    printf( "[IO] I don't understand VPF type %s\n" , token ) ;
    return FAILURE ;
  }
  printf( "[IO] We are performing %s contractions \n" , token ) ;
  return SUCCESS ;
}

// get the propagator types
static int
get_proptype( proptype *prop ,
	      const char *token ) 
{
  if( are_equal( token , "CHIRAL" ) ) {
    *prop = CHIRAL ;
  } else if( are_equal( token , "NREL" ) ) {
    *prop = NREL ;
  } else {
    printf( "[IO] I don't understand gamma basis %s\n" , 
	    token ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

///
static int
get_sourcetype( sourcetype *source ,
		const char *token ) 
{
  if( are_equal( token , "WALL" ) ) {
    *source = WALL ;
  } else if( are_equal( token , "POINT" ) ) {
    *source = POINT ;
  } else {
    printf( "[IO] I don't understand source type %s\n" , token ) ;
    return FAILURE ;
  }
  printf( "[IO] Propagators are %s sources \n" , token ) ;
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
		 const int meas_idx ) 
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
  printf( "[IO] Meson_%d :: Contracting prop %d with prop %d \n" , 
	  meas_idx , mesons -> map[0] , mesons -> map[1] ) ;

  // check for sourcetype
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_sourcetype( &( mesons -> source ) , token ) == FAILURE ) return FAILURE ;

  // check for proptype
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( mesons -> proptype1 ) , token ) == FAILURE ) return FAILURE ;
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( mesons -> proptype2 ) , token ) == FAILURE ) return FAILURE ;

  // output file
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  sprintf( mesons -> outfile , "%s" , token ) ;
  printf( "[IO] Meson_%d :: Contraction file in %s \n" , meas_idx , token ) ;

  // tell us if we get more than we expect
  if( ( token = (char*)strtok( NULL , "," ) ) != NULL ) {
    printf( "[IO] Meson_%d :: Unexpected extra contraction info %s \n" ,
	    meas_idx , token ) ;
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
			   nprops , *nmesons ) 
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

  // check for sourcetype
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_current_type( &( VPF -> current ) , token ) == FAILURE ) return FAILURE ;

  // check for proptype
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( VPF -> proptype1 ) , token ) == FAILURE ) return FAILURE ;
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( VPF -> proptype2 ) , token ) == FAILURE ) return FAILURE ;

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

  // check for proptype
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( wme -> proptype1 ) , token ) == FAILURE ) return FAILURE ;
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( wme -> proptype2 ) , token ) == FAILURE ) return FAILURE ;
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( wme -> proptype3 ) , token ) == FAILURE ) return FAILURE ;
  if( ( token = (char*)strtok( NULL , "," ) ) == NULL ) return unexpected_NULL( ) ;
  if( get_proptype( &( wme -> proptype4 ) , token ) == FAILURE ) return FAILURE ;

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
  char str[ 32 ] ;
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    sprintf( str , "DIMS%d" , mu ) ;
    const int dims_idx = tag_search( str ) ;
    if( dims_idx == FAILURE ) { 
      return tag_failure( str ) ; 
    }
    dims[ mu ] = (int)atoi( INPUT[ dims_idx ].VALUE ) ;
  }
  return SUCCESS ;
}

// get prop tags
static int
get_props( struct propfile *propfiles ,
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
      strcpy( propfiles[ *nprops ].filename , INPUT[ prop_idx ].VALUE ) ;
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
  free( inputs.prop_files ) ;
  free( inputs.VPF ) ;
  free( inputs.wme ) ;
  return ;
}

// pack the input_info struct with all our contractions
int
get_input_data( struct input_info *inputs ,
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
  get_props( inputs -> prop_files , &( inputs -> nprops ) , GLU_TRUE ) ;
  if( inputs -> nprops == 0 ) { 
    printf( "[IO] No propagator files specified \n" ) ;
  }
  inputs -> prop_files = (struct propfile*)malloc( ( inputs -> nprops ) * sizeof( struct propfile ) ) ;
  get_props( inputs -> prop_files , &( inputs -> nprops ) , GLU_FALSE ) ;

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
  if( matrix_element_contractions( inputs -> wme , &( inputs -> nVPF ) , 
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

