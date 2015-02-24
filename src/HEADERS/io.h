/**
   @file io.h
   @brief IO reader
 */
#ifndef IO_H
#define IO_H

/**
   @fn void check_checksum( FILE *fprop )
*/
int
check_checksum( FILE *fprop ) ;

/**
   @fn int read_prop( FILE *fprop, struct spinor *S , const proptype prop )
   @brief read the propagator for a timeslice
   @param fprop :: propagator file
   @param S :: spinor
   @param prop :: propagator type
   @return #SUCCESS or #FAILURE
 */
int
read_prop( FILE *fprop, 
	   struct spinor *S ,
	   const proptype prop ) ;

#endif
