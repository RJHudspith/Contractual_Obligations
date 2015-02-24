/**
   @file io.h
   @brief IO reader
 */
#ifndef IO_H
#define IO_H

/**
<<<<<<< HEAD
   @fn void check_checksum( FILE *fprop )
=======
   @fn void read_infile( const int argc , const char *argv[] )
   @brief read the input file
 */
void
read_infile( const int argc , 
	     const char *argv[] ) ;

/**
   @fn void check_checksum( FILE *fprop , long int header )
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1
   @brief have a look at the checksum
   @return #SUCCESS or #FAILURE
 */
<<<<<<< HEAD
int
check_checksum( FILE *fprop ) ;
=======
void 
check_checksum( FILE *fprop , 
		long int header ) ;
>>>>>>> 1ab5d94b3290dcf61052aaf15ee2dc2436d2d3b1

/**
   @fn int read_prop( FILE *fprop, struct spinor *S , const long int header, const int tslice )
   @brief read the propagator for a timeslice
 */
int
read_prop( FILE *fprop, 
	   struct spinor *S , 
	   const long int header, 
	   const int tslice ) ;


/**
     @fn int read_nrprop( FILE *fprop, struct spinor *S , const long int header, const int tslice )
     @brief read the NRQCD propagator for a timeslice and transform it into the chiral basis
       
  */
int
read_nrprop( FILE *fprop,
       struct spinor *S ,
       const long int header,
       const int tslice ) ;


/**
     @fn int read_chiral2nrel( FILE *fprop, struct spinor *S , const long int header, const int tslice )
     @brief read the light propagator for a timeslice and transform it into the non-rel. basis
       
  */
int
read_chiral2nrel( FILE *fprop,
       struct spinor *S ,
       const long int header,
       const int tslice ) ;




#endif
