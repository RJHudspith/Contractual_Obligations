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
   @fn void check_checksum( FILE *fprop )
   @brief have a look at the checksum
 */
void 
check_checksum( FILE *fprop ) ;

/**
   @fn int read_prop( FILE *fprop, struct spinor *S )
   @brief read the propagator for a timeslice
 */
int
read_prop( FILE *fprop, 
	   struct spinor *S ) ;


/**
     @fn int read_nrprop( FILE *fprop , struct spinor *S , const int tslice )
     @brief read the NRQCD propagator for a timeslice and transform it into the chiral basis
       
  */
int
read_nrprop( FILE *fprop,
	     struct spinor *S ) ;



#endif
