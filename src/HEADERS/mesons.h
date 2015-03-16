/**
   @file wall_mesons.h
   @brief wall meson contractions
 */

#ifndef WALL_MESONS_H
#define WALL_MESONS_H

/**
   @fn int mesons_diagonal( struct propagator prop , const char *outfile )
   @brief (flavour diagonal) meson contractions
 */
int
mesons_diagonal( struct propagator prop ,
		 const char *outfile ) ;

/**
   @fn int wall_double_mesons( struct propagator prop1 , struct propagator prop2 const char *outfile )
   @brief (flavour off diagonal) meson contractions
 */
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const char *outfile ) ;

#endif
