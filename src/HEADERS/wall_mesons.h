/**
   @file wall_mesons.h
   @brief wall meson contractions
 */

#ifndef WALL_MESONS_H
#define WALL_MESONS_H

/**
   @fn int wall_mesons( struct propagator prop , const char *outfile )
   @brief wall meson code
 */
int
wall_mesons( struct propagator prop ,
	     const char *outfile ) ;

/**
   @fn int wall_double_mesons( struct propagator prop1 , struct propagator prop2 const char *outfile )
   @brief wall meson code for two wall props
 */
int
wall_double_mesons( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const char *outfile ) ;

#endif
