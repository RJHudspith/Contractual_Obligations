/**
   @file cut_output.h
   @brief prototype functions for the outputting of momentum space data
 */

#ifndef CUT_OUTPUT_H
#define CUT_OUTPUT_H

/**
   @fn void write_projected_data( const char *filename , const int *__restrict NMOM , const double *__restrict data , const struct veclist *__restrict list , const int DIR )
   @brief write out momentum-space data
 */
void
write_momspace_data( const char *filename ,  
		     const int *__restrict NMOM ,
		     const double *__restrict data ,
		     const struct veclist *__restrict list ,
		     const int DIR ) ;

/**
   @fn void write_tmoments( const double **tcorr , const char *filename , const current_type current , const vector_axial VA )
   @brief write out the (real part of) the temporal correlators
 */
void
write_tmoments( const double **tcorr ,
		const char *filename ,
		const current_type current ,
		const vector_axial VA ) ;

#endif
