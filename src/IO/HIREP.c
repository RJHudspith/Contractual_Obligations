/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (HIREP.c) is part of GLU.

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
   @file HIREP.c
   @brief code to read and write HIREP configurations

  bizarrely their geometry is
  t,x,y,z which makes little sense to me. I imagine there is a very good reason
  for this.... Probably.
 */
#include "common.h"

#include "geometry.h"    // for changing to our geometry 
#include "GLU_bswap.h"   // byte swapping arrays
#include "plaqs_links.h" // compute the plaquette

// takes a HiRep index "i" and converts to the Nersc or my index
static size_t
translate_to_GLU_geom( const size_t i )
{
  int x[ ND ] ;
  get_mom_2piBZ_hirep2( x , i ) ;
  return gen_site( x ) ;
}

// generic (ish) reader for files in the HiRep format
int
read_gauge_field( struct site *__restrict lat , 
		  FILE *__restrict in )
{
  double Nelements = NCNC ;
#if NC==2
  Nelements = NC ;
#endif
  const size_t stride = 2 * Nelements * ND ;
  // these guys also seem to save only in double, making things easy
  double *uind = malloc( stride * sizeof( double ) ) ; 

  size_t i ;
  for( i = 0 ; i < Latt.Volume ; i++ ) {
    const size_t idx = translate_to_GLU_geom( i ) ;

    if( fread( uind , sizeof( double ) , stride , in ) != stride ) {
      fprintf( stderr , "[IO] HiRep File read error.. Leaving \n " ) ;
      free( uind ) ;
      return FAILURE ;
    }
    if ( !WORDS_BIGENDIAN ) {
      bswap_64( stride , uind ) ; 
    }

#if NC==2
    size_t mu , a = 0 ;
    // t-direction is first?
    lat[idx].O[ND-1][0] =   uind[a+0] + I * uind[a+3] ;
    lat[idx].O[ND-1][1] =  -uind[a+2] + I * uind[a+1] ;
    lat[idx].O[ND-1][2] = -conj( lat[idx].O[ND-1][1] ) ;
    lat[idx].O[ND-1][3] =  conj( lat[idx].O[ND-1][0] ) ;
    a += 2*NC ;
    // then the others
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      lat[idx].O[mu][0] =  uind[a+0] + I * uind[a+3] ;
      lat[idx].O[mu][1] =  -uind[a+2] + I * uind[a+1] ;
      lat[idx].O[mu][2] = -conj( lat[idx].O[mu][1] ) ;
      lat[idx].O[mu][3] =  conj( lat[idx].O[mu][0] ) ;
      a += 2*NC ;
    }
#else
    size_t mu , j , a = 0 ;
    // t first
    for( j = 0 ; j < NCNC ; j++ ) {
      lat[idx].O[ND-1][j] = uind[a] + I * uind[a + 1] ;
      a += 2 ;
    }
    // then the others
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[idx].O[mu][j] = uind[a] + I * uind[a + 1] ;
	a += 2 ;
      }
    }
#endif
  }
  free( uind ) ;
  return SUCCESS ; 
}
