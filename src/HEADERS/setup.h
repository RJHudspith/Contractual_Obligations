/**
   @file setup.h
   @brief prototype declarations for common setup code
 */
#ifndef CORR_SETUP_H
#define CORR_SETUP_H

/**
   @fn int compute_correlator( struct mcorr **corr , const double complex **in , const double complex **out , const struct veclist *list , const int *NMOM , const void *forward , const size_t stride1 , const size_t stride2 , const size_t tshifted )
   @brief compute the momentum-projected correlation function in @corr
   @return #SUCCESS or #FAILURE
 */
int
compute_correlator( struct mcorr **corr , 
		    const double complex **in , 
		    const double complex **out , 
		    const struct veclist *list ,
		    const int *NMOM ,
		    const void *forward , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted ) ;

/**
   @fn int free_ffts( double complex **in , double complex **out , void *forward , void *backward , const size_t flat_dirac )
   @brief free the fftw data if it has been allocated
   @return #SUCCESS or #FAILURE
 */
int
free_ffts( double complex **in , 
	   double complex **out , 
	   void *forward ,
	   void *backward ,
	   const size_t flat_dirac ) ;

/**
   @fn int init_moms( int **NMOM , int **wwNMOM , struct veclist **list , struct veclist **wwlist , const struct cut_info CUTINFO , const GLU_bool is_wall )
   @brief initialise and allocate momentum lists
   @return #SUCCESS or #FAILURE
 */
int
init_moms( int **NMOM , 
	   int **wwNMOM ,
	   struct veclist **list ,
	   struct veclist **wwlist ,
	   const struct cut_info CUTINFO , 
	   const GLU_bool is_wall ) ;

#endif
