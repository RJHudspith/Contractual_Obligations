/**
   @file correlators.h
   @brief correlator IO and storage types prototype functions
 */
#ifndef CORRELATORS_H
#define CORRELATORS_H

/**
   @fn struct mcorr **allocate_momcorrs( const int length1 , const int length2 , const int nmom )
   @brief allocate Corr[length1][length2].mom[ nmom ].C[ #L0 ] correlation function
 */
struct mcorr **
allocate_momcorrs( const int length1 , 
		   const int length2 ,
		   const int nmom ) ;

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
	      const struct mcorr **corr ) ;

/**
   @fn void debug_baryons( const char *message , const struct correlator **corr )
   @brief print to stdout some (baryon) correlator information
 */
void
debug_baryons( const char *message , 
	       const struct mcorr **corr ) ;

/**
   @fn void write_baryons( struct mcorr **Buud_corr , struct mcorr **Buuu_corr , struct mcorr **Buds_corr , const struct veclist *list , const int NMOM[ 1 ] , const GLU_bool is_wall , const char *outfile )
   @brief write out our baryonic, momentum-projected correlation functions
 */
void
write_baryons( struct mcorr **Buud_corr , 
	       struct mcorr **Buuu_corr ,
	       struct mcorr **Buds_corr ,
	       const struct veclist *list ,
	       const int NMOM[ 1 ] ,
	       const GLU_bool is_wall ,
	       const char *outfile ) ;

/**
   @fn void write_momcorr( const char *outfile , const struct mcorr **corr , const struct veclist *list , const int NSRC , const int NSNK , const int *nmom )
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
