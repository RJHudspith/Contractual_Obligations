/*
    Copyright 2013 Renwick James Hudspith

    This file (cut_routines.h) is part of GLU.

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
   @file cut_routines.h
   @brief Cutting of extreme momenta from our data
 */
#ifndef GLU_CUT_ROUTINES_H
#define GLU_CUT_ROUTINES_H

/**
  @param rats :: Anisotropy ratios
 */
extern double rats[ ND ] ;

/**
   @fn void simorb_ratios( const int DIMS , const GLU_bool configspace )
   @brief compute the ratio of the smallest to the largest lattice direction
   @param DIMS :: dimensions of the problem

   This function computes the ratios of anisotropy so that we can keep
   momenta on an equivalent standing. This is very important for the 
   non-exceptional schemes as otherwise we would not conserve the correct 
   momenta.
 **/
void
simorb_ratios( const int DIMS , 
	       const GLU_bool configspace ) ;

/**
   @fn struct veclist* compute_veclist( int *list_size , const struct cut_info CUTINFO , const int DIMS , const GLU_bool CONFIGSPACE ) ;
   @brief compute the list of momenta, or read it if possible
   @param list_size :: size of the veclist struct
   @param CUTINFO :: momentum cut information
   @param DIMS :: dimensions of the problem
   @param CONFIGSPACE :: Qsusc and Statpot have different momenta look-up
 */
struct veclist*
compute_veclist( int *list_size , 
		 const struct cut_info CUTINFO ,
		 const size_t DIMS ,
		 const GLU_bool CONFIGSPACE ) ;

/**
   @fn struct veclist_int* compute_veclist_int( int *list_size , const struct cut_info CUTINFO , const int DIMS , const GLU_bool CONFIGSPACE ) ;
   @brief compute the list of momenta, or read it if possible
   @param list_size :: size of the veclist struct
   @param CUTINFO :: momentum cut information
   @param DIMS :: dimensions of the problem
   @param CONFIGSPACE :: Qsusc and Statpot have different momenta look-up
 */
struct veclist_int*
compute_veclist_int( int *list_size , 
		     const struct cut_info CUTINFO ,
		     const size_t DIMS ,
		     const GLU_bool CONFIGSPACE ) ;

/**
   @fn struct veclist* zero_veclist( int *list_size , const int DIMS , const GLU_bool CONFIGSPACE )
   @brief returns a zeroed veclist structure
 */
struct veclist*
zero_veclist( int *list_size ,
	      const size_t DIMS ,
	      const GLU_bool CONFIGSPACE ) ;

/**
   @fn struct veclist* wall_mom_veclist( int *list_size , const double sum_mom[ ND-1 ] , const int DIMS )
   @brief returns a two-element list of momenta for the DFT
 */
struct veclist*
wall_mom_veclist( int *list_size ,
		  const double sum_mom[ ND-1 ] ,
		  const int DIMS ) ;

struct veclist*
DFT_mom_veclist( int *list_size ,
		 const struct cut_info CUTINFO ,
		 const int DIMS ) ;

#endif
