/**
   @file cut_output.h
   @brief prototype functions for the outputting of momentum space data
 */
#ifndef CUT_OUTPUT_H
#define CUT_OUTPUT_H

/**
   @fn void write_mom_veclist( FILE *Ap , const double twist[ ND ] , const int *num_mom , const struct veclist *list , const int DIR )
   @brief write out the momentum list stored in list
 */
void
write_mom_veclist( FILE *Ap ,
		   const double twist[ ND ] ,
		   const int *num_mom , 
		   const struct veclist *ist ,
		   const int DIR ) ;
/**
   @fn void write_momspace_data( const char *filename , const int *NMOM , const double *data , const struct veclist *list , const int DIR )
   @brief write out momentum-space data
 */
void
write_momspace_data( const char *filename ,  
		     const int *NMOM ,
		     const double *data ,
		     const struct veclist *list ,
		     const int DIR ) ;

#endif
