/**
   @file setup.h
   @brief prototype declarations for common setup code
 */
#ifndef CORR_SETUP_H
#define CORR_SETUP_H

/**
   @fn int compute_correlator2( struct measurements *M , const size_t stride1 , const size_t stride2 , const size_t tshifted )
   @brief compute the momentum-projected correlation function in @M.corr
 */
int
compute_correlator( struct measurements *M , 
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

/**
   @fn void free_measurements( struct measurements *M , const size_t Nprops , const size_t stride1 , const size_t stride2 , const size_t flat_dirac )
   @brief free our measurement struct and all the stuff within it
 */
void
free_measurements( struct measurements *M ,
		   const size_t Nprops ,
		   const size_t stride1 , 
		   const size_t stride2 , 
		   const size_t flat_dirac ) ;

/**
   @fn int init_measurements( struct measurements *M , const struct propagator *prop , const size_t Nprops , const struct cut_info CUTINFO , const size_t stride1 , const size_t stride2 , const size_t flat_dirac )
   @brief initialise our measurement struct
   @return #SUCCESS or #FAILURE
 */
int
init_measurements( struct measurements *M ,
		   const struct propagator *prop ,
		   const size_t Nprops ,
		   const struct cut_info CUTINFO ,
		   const size_t stride1 ,
		   const size_t stride2 ,
		   const size_t flat_dirac ) ;

#endif
