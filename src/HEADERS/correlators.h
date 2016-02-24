/**
   @file correlators.h
   @brief correlator IO and storage types prototype functions
 */
#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn struct mcorr **allocate_momcorrs( const size_t length1 , const size_t length2 , const size_t nmom )
   @brief allocate Corr[length1][length2].mom[ nmom ].C[ #LT ] correlation function
 */
struct mcorr **
allocate_momcorrs( const size_t length1 , 
		   const size_t length2 ,
		   const size_t nmom ) ;

/**
   @fn void free_momcorrs( struct mcorr **mcorr , const size_t length1 , const size_t length2 , const size_t nmom )
   @brief free the allocated mcorr struct
 */
void
free_momcorrs( struct mcorr **mcorr , 
	       const size_t length1 ,
	       const size_t length2 ,
	       const size_t nmom ) ;

/**
   @fn void debug_mesons( const char *message , const struct correlator **corr )
   @brief print to stdout some correlator information
 */
void
debug_mesons( const char *message , 
	      const struct mcorr **corr ) ;

/**
   @fn void debug_baryons( const char *message , const struct correlator **corr )
   @brief print to stdout some (baryon) correlator information
 */
void
debug_baryons( const char *message , 
	       const struct mcorr **corr ) ;

/**
   @fn void write_momcorr( const char *outfile , const struct mcorr **corr , const struct veclist *list , const size_t NSRC , const size_t NSNK , const int *nmom , const char *type )
   @brief write out the #ND-1 momentum-injected correlator
 */
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const size_t NSRC ,
	       const size_t NSNK ,
	       const int *nmom , 
	       const char *type ) ;

#endif
