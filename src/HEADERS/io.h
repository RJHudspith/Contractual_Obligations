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
   @fn int read_prop( FILE *fprop, struct spinor *S )
   @brief read the propagator for a timeslice
 */
int
read_prop( FILE *fprop, 
	   struct spinor *S ) ;


/**
     @fn int read_nrprop( FILE *fprop, struct spinor *S )
     @brief read the NRQCD propagator for a timeslice and transform it into the chiral basis
  */
int
read_nrprop( FILE *fprop,
	     struct spinor *S ) ;


/**
     @fn int read_chiral2nrel( FILE *fprop, struct spinor *S )
     @brief read the light propagator for a timeslice and transform it into the non-rel. basis
       
  */
int
read_chiral2nrel( FILE *fprop,
		  struct spinor *S ) ;




#endif
