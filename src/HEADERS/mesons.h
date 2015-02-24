/**
   @file mesons.h
   @brief prototype functions for meson computations
 */

#ifndef MESONS_H
#define MESONS_H

/**
   @fn void allocate_corrs( struct correlator **corr )
   @brief allocate NS*NS correlator matrix
 */
void
allocate_corrs( struct correlator **corr ) ;

/**
   @fn void free_corrs( struct correlator **corr )
   @brief free correlator matrix
 */
void
free_corrs( struct correlator **corr ) ;

#ifdef DEBUG
/**
   @void void debug_mesons( const char *message , const struct correlator **corr )
   @brief print to stdout some correlator information
 */
void
debug_mesons( const char *message , 
	      const struct correlator **corr ) ;
#endif

/**
   @fn int single_mesons( FILE *fprop , const proptype proptype1 )
   @brief compute mesons from a single propagator

   @return #SUCCESS or #FAILURE
 */
int
single_mesons( FILE *fprop ,
	       const proptype proptype1 ) ;

/**
   @fn int double_mesons( FILE *prop1 , const proptype proptype1 , FILE *prop2 , const proptype proptype2 )
   @brief compute mesons from two propagators

   @return #SUCCESS or #FAILURE
 */
int
double_mesons( FILE *prop1 , 
	       const proptype proptype1 ,
	       FILE *prop2 ,
	       const proptype proptype2 ) ;

#endif
