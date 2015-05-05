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
   @fn struct mcorr **allocate_momcorrs( const int length1 , const int length2 , const int nmom )
   @brief allocate Corr[length1][length2].mom[ nmom ].C[ #L0 ] correlation function
 */
struct mcorr **
allocate_momcorrs( const int length1 , 
		   const int length2 ,
		   const int nmom ) ;

/**
   @fn void free_corrs( struct correlator **corr , const int NSRC , const int NSNK )
   @brief free correlator matrix
 */
void
free_corrs( struct correlator **corr ,
	    const int NSRC ,
	    const int NSNK ) ;

/**
   @fn void free_momcorrs( struct mcorr **mcorr , const int length1 , const int length2 , const int nmom )
   @brief free the allocated mcorr struct
 */
void
free_momcorrs( struct mcorr **mcorr , 
	       const int length1 ,
	       const int length2 ,
	       const int nmom ) ;

/**
   @fn void debug_mesons( const char *message , const struct correlator **corr )
   @brief print to stdout some correlator information
 */
void
debug_mesons( const char *message , 
	      const struct correlator **corr ) ;

/**
   @fn void debug_baryons( const char *message , const struct correlator **corr )
   @brief print to stdout some (baryon) correlator information
 */
void
debug_baryons( const char *message , 
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

/**
   @fn void write_momcorr( const char *outfile , const struct mcorr **corr , const struct veclist *list , const int NSRC , const int NSNK , const int nmom )
   @brief write out the #ND-1 momentum-injected correlator
 */
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const int NSRC ,
	       const int NSNK ,
	       const int *nmom ) ;

#endif
