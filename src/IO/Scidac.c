/*
    Copyright 2013 Renwick James Hudspith

    This file (Scidac.c) is part of GLU.

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
   @file Scidac.c
   @brief Scidac and ILDG configuration header readers and writers
 */
#include "common.h"

#include "GLU_bswap.h" // for the byteswaps
#include "GLU_timer.h" // for the date
#include "XML_info.h"  // gets the important info from the xml header info

// again, define the debug if you would like to look at the header info

// < Message Begin Mask (Internal)
#define MB_MASK ((unsigned char)0x80)

// < Message End Mask (Internal)
#define ME_MASK ((unsigned char)0x40)

#define MAX_HDR64 18

#define GLU_EOF 1 // SUCCESS is 0, GLU_FAILURE is -1

// this contains all the header data
static union {
  uint64_t int64[ MAX_HDR64 ] ;
  uint32_t int32[ 2*MAX_HDR64 ] ;
  uint16_t int16[ 4*MAX_HDR64 ] ;
  unsigned char uchr[ 8*MAX_HDR64 ] ;
} header ;

// some padding
static inline int 
lime_padding( const size_t nbytes ){ return ( nbytes%8 != 0 ) ? 8 - (nbytes%8) : 0 ; }
// accessors and what have you for the header union
static inline uint64_t
header_datalength( void ) { return header.int64[1] ; }
static inline uint32_t
magic_number( void ) { return header.int32[0] ; }
static inline uint16_t
header_version( void ) { return header.int16[2] ; }
#ifdef DEBUG_ILDG
static inline unsigned char
header_mbme( void ) { return header.uchr[6] ; }
#endif

static unsigned char *lime_hdr_rec_type = &header.uchr[16] ;

/**
   @fn static int parse_SCIDAC_hdr( FILE *infile , struct head_data *HEAD_DATA , const int record )
   @brief looks through the scidac header setting up file IO stuff
 */
static int
parse_SCIDAC_hdr( FILE *infile , 
		  struct head_data *HEAD_DATA , 
		  const int record )
{
#ifdef DEBUG_ILDG
  fprintf( stdout , "\n[IO] Header for record :: %d \n\n" , record ) ;
#endif
  // assumes we have opened it fine and are at the start of the file
  const char myname[] = "GLU::parse_SCIDAC_hdr";

  // fread the whole header ...
  const int status = fread( (void*)header.int64 , sizeof( int64_t ) , MAX_HDR64 , infile ) ;
  if( status != MAX_HDR64 ) { return GLU_EOF ; }
  
  uint32_t i_magic_no[1] = { magic_number( ) } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , i_magic_no ) ; }

#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] Magic number :: %x \n" , i_magic_no[0] ) ;
#endif
  if( i_magic_no[0] != 1164413355 ) {
    fprintf( stderr , "[IO] Scideac %s wrong magic number %x vs. %x \n" , 
	     myname , i_magic_no[0] , 1164413355 ) ;
    return FAILURE ;
  }

  // now we look at the version number
  unsigned int i_version[1] = { header_version( ) } ;
  if( !WORDS_BIGENDIAN ) { bswap_32( 1 , i_version ) ; }
#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] Reading Scidac header version %x \n" , 
	   i_version[0] ) ;
#endif

  // I don't really care about the ME or MB record type
#ifdef DEBUG_ILDG
  GLU_bool i_MB, i_ME;
  if( header_mbme( ) & MB_MASK ) { 
    i_MB = GLU_TRUE ; 
  } else {
    i_MB = GLU_FALSE ;
  }
  fprintf( stdout , "[IO] %s MB flag %d \n" , myname , i_MB ) ;
  // same for ME
  if( header_mbme( ) & ME_MASK ) { 
    i_ME = GLU_TRUE ; 
  } else {
    i_ME = GLU_FALSE ;
  }
  fprintf( stdout , "[IO] %s ME flag %d \n" , myname , i_ME ) ;
#endif

  uint64_t i_data_length[1] = { header_datalength() } ;
  if( !WORDS_BIGENDIAN ) { bswap_64( 1 , i_data_length ) ; }
#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] %s DATALENGTH %llu \n" , 
	   myname , (unsigned long long)i_data_length[0] ) ;
#endif

  int padding = lime_padding( (size_t)i_data_length[0] ) ;
#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] Using %d bytes of padding for this record \n" , 
	   padding ) ;
#endif

  // again, typebuf is uninteresting ...
  unsigned char *typebuf = (unsigned char*)lime_hdr_rec_type ;
#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] %s type %s \n" , myname , typebuf ) ;
#endif

  // if we hit the binary data we exit for a binary read, could pass the length
  // of the data but the header has already given us the logical dimensions
  if( strcmp( (const char*)typebuf , "scidac-binary-data" ) == 0 ||
      strcmp( (const char*)typebuf , "ildg-binary-data" ) == 0 ) {
    // and so we exit ...
    return GLU_EOF ; // is not end of file, but not error either
  }

  // we should be able to read the remaining data
  char *io_data = malloc( ( i_data_length[0] + 1 ) * sizeof( char ) ) ;
  if( fread( io_data , sizeof( char ) , i_data_length[0] , infile ) != i_data_length[0] ) {
    fprintf( stderr , "[IO] xml information reading failed ... Leaving \n" ) ;
    free( io_data ) ;
    return FAILURE ;
  }
  io_data[ i_data_length[0] ] = '\0' ;
  // BQCD put their plaquette and cksum at the end without xml tags hmpf!
  if( strcmp( (const char*)typebuf , "bqcd-plaquette" ) == 0 ) {
    struct head_data tmp = *HEAD_DATA ;
    tmp.plaquette = atof( io_data ) ;
    *HEAD_DATA = tmp ;
  }
  // read in the possible checksum
  if( strcmp( (const char*)typebuf , "bqcd-cksum" ) == 0 ) {
    struct head_data tmp = *HEAD_DATA ;
    unsigned long int cksum ;
    sscanf( io_data , "%lu" , &cksum ) ;
    tmp.checksum = (uint32_t)cksum ;
    *HEAD_DATA = tmp ;
  }
#ifdef DEBUG_ILDG
  fprintf( stdout , "[IO] %s \n" , io_data ) ;
#endif

  parse_and_set_xml_SCIDAC( io_data , HEAD_DATA ) ;

  // free for now
  free( io_data ) ;

  // and skip the file along by "padding" amount
  if( padding > 0 ) {
    char pad_str[ padding ] ;
    if( fread( pad_str , sizeof( char ) , padding , infile ) != (size_t)padding ) {
      fprintf( stderr , "[IO] Reading of padded data failed \n" ) ;
      return FAILURE ;
    }
  }
  return SUCCESS ;
}

// scidac file has a header, plus some xml data
// all in big endian, which is nice, and is called for the ILDG too
int
get_header_data_SCIDAC( FILE *infile ,
			struct head_data *HEAD_DATA )
{
  // loop the records, zealously set max number to be 12
  int record = 0 , validity ;
  while( record < 12 ) { // hmmm
    validity = parse_SCIDAC_hdr( infile , HEAD_DATA , record++ ) ;
    if( validity == FAILURE ) { // if the reader doesn't work
      return FAILURE ;
    } else if( validity == GLU_EOF ) {
      break ;
    }
  }
#ifdef DEBUG_ILDG
  struct head_data tem = *HEAD_DATA ;
  fprintf( stdout , "[IO] ENDIAN :: %d \n" , tem.endianess ) ;
  fprintf( stdout , "[IO] PRECISION :: %d \n" , tem.precision ) ;
  fprintf( stdout , "[IO] CONFIG :: %d \n" , tem.config_type ) ;
#endif
  return SUCCESS ;
}

// clean this up
#ifdef DEBUG_ILDG
  #undef DEBUG_ILDG
#endif

// and undefs to clear up local macros
#undef A_BIG_NUMBER 
#undef MB_MASK
#undef ME_MASK
#undef MAX_HDR64
#undef GLU_EOF
