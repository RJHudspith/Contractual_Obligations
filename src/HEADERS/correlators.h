/**
   @file correlators.h
   @brief correlator IO and storage types prototype functions
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn void allocate_corrs( struct correlator **corr )
   @brief allocate \f$ NS*NS \f$ correlator matrix
 */
void
allocate_corrs( struct correlator **corr ) ;

/**
   @fn void free_corrs( struct correlator **corr )
   @brief free correlator matrix
 */
void
free_corrs( struct correlator **corr ) ;

/**
   @fn void debug_mesons( const char *message , const struct correlator **corr )
   @brief print to stdout some correlator information
 */
void
debug_mesons( const char *message , 
	      const struct correlator **corr ) ;

#endif
