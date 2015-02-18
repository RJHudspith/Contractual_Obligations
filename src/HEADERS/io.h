/**
   @file io.h
   @brief IO reader
 */
#ifndef IO_H
#define IO_H

/**
   @fn void read_infile( const int argc , const char *argv[] )
   @brief read the input file
 */
void
read_infile( const int argc , 
	     const char *argv[] ) ;

/**
   @fn void check_checksum( FILE *fprop , long int header )
   @brief have a look at the checksum
 */
void 
check_checksum( FILE *fprop , 
		long int header ) ;

/**
   @fn int read_prop( FILE *fprop, struct spinor *S , const long int header, const int tslice )
   @brief read the propagator for a timeslice
 */
int
read_prop( FILE *fprop, 
	   struct spinor *S , 
	   const long int header, 
	   const int tslice ) ;

#endif
