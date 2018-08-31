/**
   @file read_propheader.c
   @brief propagator header reader

   The way this works is the code reads a line associated with a "tag"

   e.g.
   Lattice: {stuff}

   If the tag is something we want it calls a get_* function
   these then increment whitespace until the end of the line looking for
   what we wish in the case of the lattice it will look for data of the form

   DIM_X DIM_Y DIM_Z DIM_T and compare to the global geometry

   In the case of places where a single value is expected it returns
   as soon as it has found the value

   e.g.
   Endian: Big Little Medium
   Would return Big
 */
#include "common.h"

#include <errno.h>    // for the error codes to strto*

// dirty string equivalence
static int
are_equal( const char *pch , const char *tag )
{
  char str[ strlen( tag ) + 3 ] ;
  sprintf( str , "%s\n" , tag ) ;
  return !strcmp( pch , tag ) || !strcmp( pch , str ) ;
}

// get a particular tag
static char*
get_tag( char line[ MAX_LINE_LENGTH ] )
{
  return strtok( line , " " ) ;
}

static int
skip_trailing_spaces( void )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) { }
  return SUCCESS ;
}

static int
get_double( double *result )
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    result[ N ] = strtod( token , &endptr ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] misread double %zu read %f \n" , 
	       N , result[ N ] ) ;
      return FAILURE ;
    }
    if( ++N == 1 ) return SUCCESS ;
  }
  return FAILURE ;
}

static int
get_GLU_bool( GLU_bool *result )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    if( are_equal( token , "GLU_TRUE" ) ) {
      *result = GLU_TRUE ;
      return skip_trailing_spaces( ) ;
    } else if( are_equal( token , "GLU_FALSE" ) ) {
      *result = GLU_FALSE ;
      return skip_trailing_spaces( ) ;
    } else {
      fprintf( stderr , "[IO] I don't understand bool %s, please make"
	       " GLU_FALSE or GLU_TRUE\n" , token ) ;
      return FAILURE ;
    }
  }
  return FAILURE ;
}

// get the initial source position
// this is really gross as strcmp doesn't play nice with the hyphen or a
// newline so I split the string again and look for anti
// also there appears to be different capitalisations in the props!
static int
get_propbounds( boundaries *bound ) 
{
  char *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    if( are_equal( token , "periodic" ) || are_equal( token , "Periodic" ) ) {
      bound[ N ] = PERIODIC ;
    } else if( are_equal( token , "PplusA" ) ) {
      bound[ N ] = PPLUSA ;
    } else if( are_equal( token , "PminusA" ) ) {
      bound[ N ] = PMINUSA ;
    } else if( are_equal( token , "PmulA" ) ) {
      bound[ N ] = PMULA ;
    } else {
      char *tok2 = strtok( token , "-" ) ;
      if( are_equal( "anti" , tok2 ) || are_equal( "Anti" , tok2 )) {
	bound[ N ] = ANTIPERIODIC ;
      } else {
	fprintf( stderr , "[IO] propheader I don't understand boundary %s\n" ,
		 token ) ;
	return FAILURE ;
      }
    }
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

// get the propagator dims
static int
get_propdims( size_t *dims ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    dims[ N ] = (size_t)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] propheader dims[%zu] misread %zu \n" , 
	       N , dims[ N ] ) ;
      return FAILURE ;
    }
    // check against global lattice data struct
    if( dims[ N ] != Latt.dims[ N ] ) {
      fprintf( stderr , "[IO] propheader and global lattice dims mismatch!\n"
	       "[IO] %zu vs %zu ( index %zu ) \n" , dims[ N ] ,
	       Latt.dims[ N ] , N ) ;
      return FAILURE ;
    }
    // we have successfully read in all the dimension info we need
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

// get the propagator endianness
static int
get_propendian( endianness *endian )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Big" ) ) {
      *endian = BIGENDIAN ;
      return SUCCESS ;
    }
    if( are_equal( token , "Little" ) ) {
      *endian = LILENDIAN ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// get the propagator precision
static int
get_propprec( fp_precision *precision )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Single" ) ) {
      *precision = SINGLE ;
      return SUCCESS ;
    }
    if( are_equal( token , "Double" ) ) {
      *precision = DOUBLE ;
      return SUCCESS ;
    }
    //
  }
  return FAILURE ;
}

// get the initial source position
static int
get_propsmear( smearing *smear )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Quark" ) ) {
      fprintf( stdout , "[IO] Smeared quark source\n" ) ;
      *smear = QUARK ;
      return SUCCESS ;
    } else if( are_equal( token , "Gauge" ) ) {
      *smear = GAUGE ;
      return SUCCESS ;
    } else {
      *smear = NOSMEAR ;
      return SUCCESS ;
    }
  }
  return FAILURE ;
}

// get the initial source position
static int
get_propsrc( size_t *origin ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    origin[ N ] = (size_t)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE || 
	origin[ N ] > Latt.dims[ N ] ) {
      fprintf( stderr , "[IO] propheader SrcPos:[%zu] misread %zu \n" , 
	       N , origin[ N ] ) ;
      return FAILURE ;
    }
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

// get the propagator source type
static int
get_propsource( sourcetype *source )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // set the endian flag
    if( are_equal( token , "Point" ) || 
	are_equal( token , "Pt" ) || 
	are_equal( token , "pt" ) ) {
      *source = POINT ;
      return SUCCESS ;
    }
    // wall can be a bunch of things I think
    if( are_equal( token , "GFWall" ) ||
	are_equal( token , "Wall" ) ||
	are_equal( token , "wall" ) || 
	are_equal( token , "CGwall" ) ) {
      *source = WALL ;
      return SUCCESS ;
    }
    if( are_equal( token , "Z2Wall" ) ||
	are_equal( token , "Z2_Wall" ) ||
	are_equal( token , "Z2_WALL" ) ) {
      *source = Z2_WALL ;
      return SUCCESS ;
    }
  }
  return FAILURE ;
}

// get the propagator basis type
static int
get_proptype( proptype *basis )
{
  char *token ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    // nonrelativistic basis
    if( are_equal( token , "Nrel_fwd" ) ) {
      *basis = NREL_FWD ;
      return skip_trailing_spaces( ) ;
    } else if( are_equal( token , "Nrel_bwd" ) ) {
      *basis = NREL_BWD ;
      return skip_trailing_spaces( ) ;
      // on the fly computed NRQCD tag
    } else if( are_equal( token , "Nrel_CORR" ) ) {
      *basis = NREL_CORR ;
      return skip_trailing_spaces( ) ;
    } else if( are_equal( token , "Chiral" ) ) {
      *basis = CHIRAL ;
      return skip_trailing_spaces( ) ;
    }
    //
  }
  return FAILURE ;
}

// get propagator twists
static int
get_proptwists( double *twists ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    twists[ N ] = strtod( token , &endptr ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] propheader twists[%zu] misread %f \n" , 
	       N , twists[ N ] ) ;
      return FAILURE ;
    }
    // we have successfully read in all the twist info we need
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

// get propagator twists
static int
get_propmom_source( double *mom_source ) 
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    mom_source[ N ] = strtod( token , &endptr ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] propheader mom_source[%zu] misread %f \n" , 
	       N , mom_source[ N ] ) ;
      return FAILURE ;
    }
    // we have successfully read in all the twist info we need
    if( ++N == ND ) return SUCCESS ;
  }
  return FAILURE ;
}

static int
get_size_t( size_t *result )
{
  char *endptr , *token ;
  size_t N = 0 ;
  errno = 0 ;
  while( ( token = strtok( NULL , " " ) ) != NULL ) {
    result[ N ] = (size_t)strtol( token , &endptr , 10 ) ;
    if( token == endptr || errno == ERANGE ) {
      fprintf( stderr , "[IO] misread size_t %zu read %zu \n" , 
	       N , result[ N ] ) ;
      return FAILURE ;
    }
    if( ++N == 1 ) return SUCCESS ;
  }
  return FAILURE ;
}

// little sugar
static int
tagfailure( const char *message , 
	    const char *line )
{
  fprintf( stderr , "[IO] propheader failure while searching for tag %s \n"
	   "[IO] line in file is %s \n" , message , line ) ;
  return FAILURE ;
}

// tells us if we missed a record we were looking for
static int
nonexistent_record( const char *message )
{
  fprintf( stderr , "[IO] I cannot find record for tag %s \n" , message ) ;
  return FAILURE ;
}

static void
summarize_NRQCD_params( struct NRQCD_params NRQCD )
{
  fprintf( stdout , "\n[IO] NRQCD bare mass %f\n" , NRQCD.M_0 ) ;
  fprintf( stdout , "[IO] NRQCD tadpole %f\n" , NRQCD.U0 ) ;
  fprintf( stdout , "[IO] NRQCD %zu applications of hamiltonian\n" , NRQCD.N ) ;
  if( NRQCD.backward == GLU_TRUE ) {
    fprintf( stdout , "[IO] NRQCD backward propagator\n" ) ;
  } else {
    fprintf( stdout , "[IO] NRQCD forward propagator\n" ) ;
  }
#ifdef NRQCD_NONSYM
  fprintf( stdout , "[IO] NRQCD single application of spin-dependent part\n" ) ;
#else
  fprintf( stdout , "[IO] NRQCD symmetric spin-dependent application\n" ) ;
#endif
#ifdef LEGACY_NRQCD_COMPARE
  fprintf( stdout , "[IO] NRQCD *NOT* subtracting tadpole terms from"
	   " improved derivatives\n" ) ;
#else
  fprintf( stdout , "[IO] NRQCD subtracting tadpole terms from"
	   " improved derivatives\n" ) ;
#endif
  fprintf( stdout , "[IO] NRQCD coefficient C0 %f\n" , NRQCD.C0 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C1 %f\n" , NRQCD.C1 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C2 %f\n" , NRQCD.C2 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C3 %f\n" , NRQCD.C3 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C4 %f\n" , NRQCD.C4 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C5 %f\n" , NRQCD.C5 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C6 %f\n" , NRQCD.C6 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C7 %f\n" , NRQCD.C7 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C8 %f\n" , NRQCD.C8 ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C9EB %f\n" , NRQCD.C9EB ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C10EB %f\n" , NRQCD.C10EB ) ;
  fprintf( stdout , "[IO] NRQCD coefficient C11 %f\n" , NRQCD.C11 ) ;
  return ;
}

static void
summarize_prop_source( const struct propagator prop )
{
  size_t i ;
  fprintf( stdout , "[IO] origin:" ) ;
  for( i = 0 ; i < ND ;  i++ ) {
    fprintf( stdout , " %zu" , prop.origin[i] ) ;
  }
  fprintf( stdout , "\n" ) ;
  
  fprintf( stdout , "[IO] twist:" ) ;
  for( i = 0 ; i < ND ;  i++ ) {
    fprintf( stdout , " %g" , prop.twist[i] ) ;
  }
  fprintf( stdout , "\n" ) ;

  fprintf( stdout , "[IO] mom source:" ) ;
  for( i = 0 ; i < ND ;  i++ ) {
    fprintf( stdout , " %g" , prop.mom_source[i] ) ;
  }
  fprintf( stdout , "\n" ) ;  

  if( prop.Source.smear == QUARK ) {
    fprintf( stdout , "[IO] prop has smeared quark source\n" ) ;
    fprintf( stdout , "[IO] prop smearing with %zu iterations\n" ,
	     prop.Source.Nsmear ) ;
    fprintf( stdout , "[IO] prop using smearing alpha %f\n" ,
	     prop.Source.smalpha ) ;
  }

  switch( prop.Source.type ) {
  case POINT :
    fprintf( stdout , "[IO] propagator is a POINT source\n" ) ;
    break ;
  case WALL :
    fprintf( stdout , "[IO] propagator is a WALL/Box source\n" ) ;
    fprintf( stdout , "[IO] Wall box size is %zu\n" , prop.Source.boxsize ) ;
    break ;
  case Z2_WALL :
    fprintf( stdout , "[IO] propagator is a Z2_WALL source\n" ) ;
    fprintf( stdout , "[IO] propagator has spacing %zu\n" ,
	     prop.Source.Z2_spacing ) ;
    break ;
  }
  
  return ;
}

//
int
read_propheader( struct propagator *prop )
{
  // dimensions in the file
  size_t dims[ ND ] ;
  const int MAX_HEADER_LINES = 64 ;
  char line[ MAX_LINE_LENGTH ] ;

  // error flags
  int n = 0 , dimsflag = 0 , originflag = 0 ;
  int endflag = 0 , srcflag = 0 , precflag = 0 ;
  int basisflag = 0 , boundsflag = 0 , twistsflag = 0 ;
  int momsourceflag = 0 , smearingflag = 0 ;
  
  // set this to NULL
  prop -> H = NULL ;
  
  // initialise NRQCD parameters regardless of if we use them
  prop -> NRQCD.C0    = 0.0 ; prop -> NRQCD.C1   = 0.0 ;
  prop -> NRQCD.C2    = 0.0 ; prop -> NRQCD.C3   = 0.0 ;
  prop -> NRQCD.C4    = 0.0 ; prop -> NRQCD.C5   = 0.0 ;
  prop -> NRQCD.C6    = 0.0 ; prop -> NRQCD.C7   = 0.0 ;
  prop -> NRQCD.C8    = 0.0 ; prop -> NRQCD.C9EB = 0.0 ; 
  prop -> NRQCD.C10EB = 0.0 ; prop -> NRQCD.C11  = 0.0 ;
  prop -> NRQCD.M_0   = 1.0 ; prop -> NRQCD.U0   = 1.0 ;
  prop -> NRQCD.N     = 0   ; prop -> NRQCD.backward = GLU_FALSE ;

  // some defaults for the smearing and Z2 stuff

  // this is an LCU wall source, logically boxsize = 1 is a point source
  prop -> Source.boxsize = LCU ; 
  prop -> Source.Nsmear = 0 ;
  prop -> Source.smalpha = 0.15 ;
  prop -> Source.Z2_spacing = 1 ;

  // initialise these to zero
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    prop -> twist[ mu ] = prop -> mom_source[ mu ] = 0.0 ;
  }
  
  // loop up to some large number of possible tags
  while( n < MAX_HEADER_LINES ) {

    // read a line
    if( fgets( line , MAX_LINE_LENGTH , prop -> file ) == NULL ) {
      fprintf( stderr , "[IO] prop Header reading failed ... Leaving \n" ) ;
      return FAILURE ;
    }
    char *tag = get_tag( line ) ;

    // tag search for Lattice
    if( are_equal( tag , "Lattice:" ) ) {
      if( get_propdims( dims ) == FAILURE ) {
	return tagfailure( "Lattice:" , line ) ;
      }
      dimsflag ++ ; 
    }
    // tag search for source origin
    if( are_equal( tag , "SrcPos:" ) ) {
      if( get_propsrc( prop -> origin ) == FAILURE ) {
	return tagfailure( "SrcPos:" , line ) ;
      }
      originflag ++ ; 
    }
    // incorporate these when the header catches up
    if( are_equal( tag , "Endian:" ) ) {
      if( get_propendian( &( prop -> endian ) ) == FAILURE ) {
	return tagfailure( "Endian:" , line ) ;
      }
      endflag ++ ;
    }
    // precision
    if( are_equal( tag , "Precision:" ) ) {
      if( get_propprec( &( prop -> precision ) ) == FAILURE ) {
	return tagfailure( "Precision:" , line ) ;
      }
      precflag ++ ;
    }
    // source
    if( are_equal( tag , "Source:" ) ) {
      if( get_propsource( &( prop -> Source.type ) ) == FAILURE ) {
	return tagfailure( "Source:" , line ) ;
      }
      srcflag ++ ;
    }
    // basis
    if( are_equal( tag , "Basis:" ) ) {
      if( get_proptype( &( prop -> basis ) ) == FAILURE ) {
	return tagfailure( "Basis:" , line ) ;
      }
      basisflag ++ ;
    }
    // propagator boundaries
    if( are_equal( tag , "Boundaries:" ) || are_equal( tag , "Boundary:" )) {
      if( get_propbounds( prop -> bound ) == FAILURE ) {
	return tagfailure( "Boundaries:" , line ) ;
      }
      boundsflag++ ;
    }
    // propagator boundaries
    if( are_equal( tag , "Twists:" ) ) {
      if( get_proptwists( prop -> twist ) == FAILURE ) {
	return tagfailure( "Twists:" , line ) ;
      }
      twistsflag++ ;
    }
    // propagator mom source
    if( are_equal( tag , "Mom_Source:" ) ) {
      if( get_propmom_source( prop -> mom_source ) == FAILURE ) {
	//return tagfailure( "Mom_Source:" , line ) ;
	fprintf( stderr , "[IO] tag Mom_source: not found in prop header. Soon this will become a sinalling error - J !\n" ) ;
      }
      momsourceflag++ ;
    }
    // propagator source Smearing
    if( are_equal( tag , "Smearing:" ) ) {
      if( get_propsmear( &prop -> Source.smear ) == FAILURE ) {
	//return tagfailure( "Smearing:" , line ) ;
	fprintf( stderr , "[IO] tag Smearing: not found in prop header. Soon this will become a sinalling error - J !\n" ) ;
      }
      smearingflag++ ;
    }
    
    // look for some NRQCD parameters
    if( are_equal( tag , "NRQCD_C0" ) ) get_double( &prop -> NRQCD.C0 ) ;
    if( are_equal( tag , "NRQCD_C1" ) ) get_double( &prop -> NRQCD.C1 ) ;
    if( are_equal( tag , "NRQCD_C2" ) ) get_double( &prop -> NRQCD.C2 ) ;
    if( are_equal( tag , "NRQCD_C3" ) ) get_double( &prop -> NRQCD.C3 ) ;
    if( are_equal( tag , "NRQCD_C4" ) ) get_double( &prop -> NRQCD.C4 ) ;
    if( are_equal( tag , "NRQCD_C5" ) ) get_double( &prop -> NRQCD.C5 ) ;
    if( are_equal( tag , "NRQCD_C6" ) ) get_double( &prop -> NRQCD.C6 ) ;
    if( are_equal( tag , "NRQCD_C7" ) ) get_double( &prop -> NRQCD.C7 ) ;
    if( are_equal( tag , "NRQCD_C8" ) ) get_double( &prop -> NRQCD.C8 ) ;
    if( are_equal( tag , "NRQCD_C9EB" ) ) get_double( &prop -> NRQCD.C9EB ) ;
    if( are_equal( tag , "NRQCD_C10EB" ) ) get_double( &prop -> NRQCD.C10EB ) ;
    if( are_equal( tag , "NRQCD_C11" ) ) get_double( &prop -> NRQCD.C11 ) ;
    if( are_equal( tag , "NRQCD_U0" ) ) get_double( &prop -> NRQCD.U0 ) ;
    if( are_equal( tag , "NRQCD_M_0" ) ) get_double( &prop -> NRQCD.M_0 ) ;
    if( are_equal( tag , "NRQCD_N" ) ) get_size_t( &prop -> NRQCD.N ) ;
    if( are_equal( tag , "NRQCD_backward" ) ) get_GLU_bool( &prop -> NRQCD.backward ) ;

    // NRQCD sources
    if( are_equal( tag , "Boxsize:" ) ) get_size_t( &prop -> Source.boxsize ) ;
    if( are_equal( tag , "Nsmear:" ) ) get_size_t( &prop -> Source.Nsmear ) ;
    if( are_equal( tag , "Smalpha:" ) ) get_double( &prop -> Source.smalpha ) ;
    if( are_equal( tag , "Z2_spacing:" ) ) get_size_t( &prop -> Source.Z2_spacing ) ;
    
    // break when we hit the desired end_header
    if( are_equal( line , "<end_header>\n" ) ) {
      break ;
    }
    n++ ;
  }

  // check the flags have been hit
  if( dimsflag == 0 ) return nonexistent_record( "Lattice:" ) ; 
  if( originflag == 0 ) return nonexistent_record( "SrcPos:" ) ; 
  if( endflag == 0 ) return nonexistent_record( "Endian:" ) ; 
  if( precflag == 0 ) return nonexistent_record( "Precision:" ) ;
  if( srcflag == 0 ) return nonexistent_record( "Source:" ) ;
  if( basisflag == 0 ) return nonexistent_record( "Basis:" ) ;
  if( boundsflag == 0 ) return nonexistent_record( "Boundaries: or Boundary:" ) ;
  // I will not complain about twists not being read
  if( n == MAX_HEADER_LINES ) return nonexistent_record( "<end header>" ) ;

  if( prop -> basis == NREL_CORR ) {
    summarize_NRQCD_params( prop -> NRQCD ) ;
  }

  summarize_prop_source( *prop ) ;

  // sanity check z2_spacing
  if( prop -> Source.type == Z2_WALL ) {
    if( prop -> Source.Z2_spacing == 0 ) {
      fprintf( stderr , "[IO] prop header non sensical Z2_spacing of 0\n"
	       "[IO] setting to 1, i.e. all points are Z2xZ2\n" ) ;
      prop -> Source.Z2_spacing = 1 ;
    }
    GLU_bool bad_spacing = GLU_FALSE ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      if( Latt.dims[mu]%prop -> Source.Z2_spacing != 0 ) {
	fprintf( stderr , "[IO] Z2_spacing %zu not a factor L_%zu -> %zu!\n" ,
		 prop -> Source.Z2_spacing , mu , Latt.dims[mu] ) ;
	bad_spacing = GLU_TRUE ;
      }
      if( bad_spacing == GLU_TRUE ) return FAILURE ;
    }
  }

  // Randy's NRQCD code counts from 1 instead of zero, shift to c-counting
  // instead of Fortran counting ->  I hate this so much
  if( prop -> basis == NREL_FWD || prop -> basis == NREL_BWD ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      prop -> origin[ mu ] -= 1 ;
    }
  }
  
  // if everything is kosher we leave successfully
  return SUCCESS ;
}

// open the files and parse the header
int
read_propheaders( struct propagator *prop ,
		  const size_t nprops )
{
  size_t i ;
  for( i = 0 ; i < nprops ; i++ ) {
    // read and check 'em
    if( read_propheader( &prop[ i ] ) == FAILURE ) {
      return FAILURE ;
    }
  }
  return SUCCESS ;
}
