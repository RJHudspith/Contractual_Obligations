/**
   @file wall_mesons.h
   @brief wall meson contractions
 */

#ifndef WALL_MESONS_H
#define WALL_MESONS_H

/**
   @fn int wall_mesons( FILE *prop1 , const int header )
   @brief singlet wall meson code
 */
int
wall_mesons( FILE *prop1 , 
	     const int header ) ;

/**
   @fn int wall_double_mesons( FILE *prop1 , FILE *prop2 ,const int header )
   @brief singlet wall meson code
 */
int
wall_double_mesons( FILE *prop1 , 
		    FILE *prop2 , 
		    const int header ) ;

#endif
