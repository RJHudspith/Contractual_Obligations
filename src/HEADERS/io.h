/**
   @file io.h
   @brief IO prototype functions
 */
#ifndef IO_H
#define IO_H

/**
   @fn void check_checksum( FILE *fprop )
   @brief check the checksum at the end of the file
   @return #SUCCESS or #FAILURE
*/
int
check_checksum( FILE *fprop ) ;

/**
   @fn int read_ahead( struct propagator *prop , struct spinor *S , int *error_code , const size_t Nprops )
   @brief read the timeslice above
   @warning should be called in an OMP parallel region
 */
int
read_ahead( struct propagator *prop ,
	    struct spinor **S , 
	    int *error_code ,
	    const size_t Nprops ,
	    const size_t t ) ;

/**
   @fn int read_prop( struct propagator prop , struct spinor *S , const size_t t )
   @brief read the propagator for a timeslice
   @param prop :: propagator file
   @param S :: spinor
   @param t :: time index
   @return #SUCCESS or #FAILURE
 */
int
read_prop( struct propagator prop ,
	   struct spinor *S ,
	   const size_t t ) ;

#endif
