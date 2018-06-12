/**
   @file correlators.h
   @brief correlator computation, IO, and storage types prototype functions
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
   @fn int compute_correlator( struct measurements *M , const size_t stride1 , const size_t stride2 , const size_t tshifted )
   @brief compute the momentum-projected correlation function in @M.corr
 */
int
compute_correlator( struct measurements *M , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted ) ;

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
