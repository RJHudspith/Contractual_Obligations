/**
   @file wall_mesons.h
   @brief wall meson contractions
 */

#ifndef WALL_MESONS_H
#define WALL_MESONS_H

/**
   @fn int wall_mesons( FILE *prop1 , const proptype proptype1 , const char *outfile )
   @brief wall meson code
 */
int
wall_mesons( FILE *prop1 ,
	     const proptype proptype1 ,
	     const char *outfile ) ;

/**
   @fn int wall_double_mesons( FILE *prop1 , const proptype proptype1 , FILE *prop2 , const proptype proptype2 , const char *outfile )
   @brief wall meson code for two wall props
 */
int
wall_double_mesons( FILE *prop1 , 
		    const proptype proptype1 ,
		    FILE *prop2 ,
		    const proptype proptype2 ,
		    const char *outfile ) ;

#endif
