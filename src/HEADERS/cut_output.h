/**
   @file cut_output.h
   @brief prototype functions for the outputting of momentum space data
 */

#ifndef CUT_OUTPUT_H
#define CUT_OUTPUT_H

/**
   @fn void write_mom_veclist( FILE *__restrict Ap , const int *__restrict num_mom , const struct veclist *__restrict list , const int DIR )
   @brief write out the momentum list stored in list
 */
void
write_mom_veclist( FILE *__restrict Ap , 
		   const int *__restrict num_mom , 
		   const struct veclist *__restrict list ,
		   const int DIR ) ;
/**
   @fn void write_momspace_data( const char *filename , const int *__restrict NMOM , const double *__restrict data , const struct veclist *__restrict list , const int DIR )
   @brief write out momentum-space data
 */
void
write_momspace_data( const char *filename ,  
		     const int *__restrict NMOM ,
		     const double *__restrict data ,
		     const struct veclist *__restrict list ,
		     const int DIR ) ;

#endif
