/*
    Copyright 2013 Renwick James Hudspith

    This file (XML_info.c) is part of GLU.

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
   @file XML_info.c
   @brief dirty parser for the xml information in Scidac and ILDG config files
 */
#include "common.h"

#include <string.h>

// this prints to stdout all of the header information ...
//#define DEBUG_ILDG

// slightly clearer (for me) string equality function
static int
are_equal( const char *pch , const char *tag )
{
  return !strcmp( pch , tag ) ;
}

/**
  These functions came from ahmidas, they claim to have had trouble parsing ILDG
  xml data due to some poor choices by ETMC, just to be safe I include them as
  my code for understanding the xml is similar to theirs
 */
inline static char *realFront(char *string) { return string + strspn(string, " \t"); }

// returns 0 if there is nothing else returns the value we
// are looking for
static int
get_int_tag( char *pch , const char *tag )
{
  const int length = strlen( tag ) ;
  char tmp[ length + 2 ] ;
  sprintf( tmp , "/%s" , tag ) ;
  // loop tokens
  int value = 0 ;
  while( ( pch = strtok( 0 , "<>" ) ) ) {
    if( !strncmp( pch , tmp , length+1 ) )
      break ;
    value = atoi( realFront( pch ) ) ;
  }
  return value ;
}

// the result is in pch
static int
get_prec_tag( char *pch , const char *tag )
{
  const int length = strlen( tag ) ;
  char tmp[ length + 2 ] ;
  sprintf( tmp , "/%s" , tag ) ;
  char prec[ 8 ] = "00" ;
  // loop tokens
  while( ( pch = strtok( 0 , "<>" ) ) ) {
    if( !strncmp( pch , tmp , length+1 ) ) break ;
    sprintf( prec , "%s" , pch ) ;
  }
  if( !strcmp( prec , "F" ) ) return FLOAT_PREC ;
  if( !strcmp( prec , "D" ) ) return DOUBLE_PREC ;
  if( atoi( realFront( prec ) ) == 32 ) return FLOAT_PREC ;
  if( atoi( realFront( prec ) ) == 64 ) return DOUBLE_PREC ;
  return FAILURE ;
}

// parse the information from the c-str of xml info
int
parse_and_set_xml_SCIDAC( char *xml_info ,
			  struct head_data *HEAD_DATA ) 
{
  // scidac is always big endian, which is nice
  HEAD_DATA -> endianess = B_ENDIAN ;

  // We use the C tokenize capabilities to parse this string
  char *pch = strtok( xml_info , "<>" ) ;

  if (strncmp( pch , "?xml" , 4 ) ) {
    #ifdef DEBUG_ILDG
    printf( "[IO] No xml data found in this record ... \n" ) ;
    #endif
    return 0 ;
  }
  int dimensions = 0 ;

  // set up the search for extra dimensions, extra space for terminating null
  const char open[4][3] = { "lx" , "ly" , "lz" , "lt" } ;
  char search[ND][3] ;
  int mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    if( mu < 4 ) {
      sprintf( search[mu] , "%s" , open[mu] ) ;
    } else {
      sprintf( search[mu] , "l%d" , mu ) ;
    }
  }
  sprintf( search[ND-1] , "%s" , open[3] ) ;

  // We've removed the XML header, now we can set up a state machine to parse the file
  while ( ( pch = strtok( 0 , "<>" ) ) ) {
    while ( ( pch = strtok( 0 , "<>" ) ) ) {
      
      // break if we are at the end of a scidac file or record
      if ( are_equal( pch , "/scidacRecord" ) ) break ;
      if ( are_equal( pch , "/scidacFile" ) ) break ;
      if ( are_equal( pch , "/ildgFormat" ) ) break ;

      // get number of spacetime dimensions
      if( are_equal( pch , "spacetime" ) ) {
	if( ( dimensions = get_int_tag( pch , "spacetime" ) ) != 0 ) {
	  if( dimensions != ND ) {
	    fprintf( stderr , "[IO] ND mismatch compiled %d vs. read %d \n" ,
		     ND , dimensions ) ;
	    return FAILURE ;
	  }
	}
	#ifdef DEBUG_ILDG
	fprintf( stdout , "[IO] spacetime :: %d \n" , dimensions ) ;
	#endif
      }

      // have a look at the colors
      if( are_equal( pch , "colors" ) ) {
	if( ( dimensions = get_int_tag( pch , "colors" ) ) != 0 ) {
	  if( dimensions != NC ) {
	    fprintf( stderr , "[IO] NC mismatch compiled %d vs. read %d \n" ,
		     NC , dimensions ) ;
	    return FAILURE ;
	  }
	}
	#ifdef DEBUG_ILDG
	fprintf( stdout , "[IO] colors :: %d \n" , dimensions ) ;
	#endif
      }

      // have a look at the storage type size
      if( are_equal( pch , "typesize" ) ) {
	if( ( dimensions = get_int_tag( pch , "typesize" ) ) != 0 ) {
	  HEAD_DATA -> config_type = OUTPUT_NCxNC ; // scidac uses this binary file output
	  if( dimensions == 2 * ND * NCNC ) {
	    HEAD_DATA -> precision = FLOAT_PREC ;
	  } else if( dimensions == ( 4 * ND * NCNC ) ) {
	    HEAD_DATA -> precision = DOUBLE_PREC ;
	  } else {
	    fprintf( stderr , "[IO] Storage type not recognised \n" ) ;
	    return FAILURE ;
	  }
	}
	#ifdef DEBUG_ILDG
	fprintf( stdout , "[IO] typesize :: %d \n" , dimensions ) ;
	#endif
      }

      // get the floating point output precision ...
      if( are_equal( pch , "precision" ) ) {
	if( ( dimensions = get_prec_tag( pch , "precision" ) ) != -1 ) {
	  if( dimensions == FLOAT_PREC ) {
	    HEAD_DATA -> precision = FLOAT_PREC ;
	  } else if( dimensions == DOUBLE_PREC ) {
	    HEAD_DATA -> precision = DOUBLE_PREC ;
	  } else {
	    fprintf( stderr , "[IO] Precision %d not understood \n" , dimensions ) ;
	    return FAILURE ;
	  }
	}
      }

      // grok the lattice dimensions
      if ( are_equal( pch, "dims" ) ) {
	while ( ( pch = strtok( 0 , "<>" ) ) ) {
	  // break up the dimensions ...
	  char *token = strtok( pch , " " ) ;
	  if( token == NULL ) return FAILURE ;
	  int idx = 0 ;
	  char *pEnd ;
	  Latt.dims[ idx++ ] = strtol( token , &pEnd , 10 ) ;
	  #ifdef DEBUG_ILDG
	  fprintf( stdout , "[IO] DIMS_%d \n" , Latt.dims[idx]) ;
	  #endif
	  while( ( token = strtok( NULL , " " ) ) != NULL ) {
	    if(  ( are_equal( token , "/dims" ) ) ) break ;
	    Latt.dims[ idx ] = strtol( token , &pEnd , 10 ) ;
	    #ifdef DEBUG_ILDG
	    fprintf( stdout , "[IO] DIMS_%d \n" , Latt.dims[idx]) ;
	    #endif
	    idx++ ;
	  } 
	  // and break the xml
	  break;
	}
	continue ;
      }

      // I just take the first checksum ...
      if( are_equal( pch , "suma" ) ) {
	while( ( pch = strtok( 0 , "<>" ) ) ) {
	  if( are_equal( pch , "/suma" ) ) break ;
	  sscanf( pch , "%x" , &(HEAD_DATA -> checksum) ) ;
	  #ifdef DEBUG_ILDG
	  fprintf( stdout , "[IO] checksum %x \n" , HEAD_DATA -> checksum ) ;
	  #endif
	}
	continue ;
      }

      ////////////////////////  ILDG SPECIFIC TAGS ////////////////////
      // If I were to make up a format it would be exactly the       //
      // same as another one apart from tiny, pointless, differences //

      // geometry, not sure about making this ND-generic much prefer the
      // scidac "dims" array
      int length = 0 ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( are_equal( pch , search[mu] ) ) {
	  if( ( length = get_int_tag( pch , search[mu] ) ) != 0 ) {
	    Latt.dims[ mu ] = length ;
            #ifdef DEBUG_ILDG
	    fprintf( stdout , "[IO] ILDG Lx :: %d \n" , length ) ;
            #endif
	  }
	  continue ;
	}
      }

      // for getting the field ... su3gauge is actually NCxNC by the look of it
      // loop tokens
      if( are_equal( pch , "field" ) ) {
	while( ( pch = strtok( 0 , "<>" ) ) ) {
	  if( are_equal( pch , "/field" ) ) break ;
	  // sometimes there is whitespace around these for some unknown reason
	  char *token = realFront( pch ) ;
	  char compare[10] ;
	  sprintf( compare , "su%dgauge" , NC ) ;
	  if( !strncmp( compare , token , 8 ) ) {
	    HEAD_DATA -> config_type = OUTPUT_NCxNC ;
	    #ifdef DEBUG_ILDG
	    fprintf( stdout , "[ILDG] configuration type %d \n" , 
		     HEAD_DATA -> config_type ) ; 
	    #endif
	  } else {
	    fprintf( stderr , "[ILDG] Expected %s, got \"%s\" \n" , 
		     compare , token ) ;
	    return FAILURE ;
	  }
	}
	continue ;
      }
    }
    break;
  }
  return SUCCESS ;
}
// YUCK

#ifdef DEBUG_ILDG
  #undef DEBUG_ILDG
#endif
