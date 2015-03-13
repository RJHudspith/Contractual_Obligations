/**
   @file correlators.h
   @brief correlator IO and storage types prototype functions
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn struct correlator ** allocate_corrs( const int NSRC , const int NSNK )
   @brief allocate NSRC * NSNK correlator matrix C(t)
 */
struct correlator **
allocate_corrs( const int NSRC , 
		const int NSNK ) ;

/**
   @fn void free_corrs( struct correlator **corr )
   @brief free correlator matrix
 */
void
free_corrs( struct correlator **corr ,
	    const int NSRC ,
	    const int NSNK ) ;

/**
   @fn void debug_mesons( const char *message , const struct correlator **corr )
   @brief print to stdout some correlator information
 */
void
debug_mesons( const char *message , 
	      const struct correlator **corr ) ;

/**
   @fn void write_correlators( const char *outfile , const struct correlator **corr )
   @brief write the full correlation matrix out in binary to outfile
 */
void
write_correlators( const char *outfile ,
		   const struct correlator **corr ,
		   const int NSRC ,
		   const int NSNK ) ;

#endif
