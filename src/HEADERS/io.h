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
   @fn int read_prop( FILE *fprop, struct spinor *S )
   @brief read the propagator for a timeslice
   @param prop :: propagator file
   @param S :: spinor
   @return #SUCCESS or #FAILURE
 */
int
read_prop( struct propagator prop ,
	   struct spinor *S ) ;

#endif
