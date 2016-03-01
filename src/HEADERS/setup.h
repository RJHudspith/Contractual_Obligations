/**
   @file setup.h
   @brief prototype declarations for common setup code
 */
#ifndef CORR_SETUP_H
#define CORR_SETUP_H

/**
   @fn void copy_props( struct measurements *M , const size_t Nprops )
   @brief set the pointers to the next timeslice's memory address
 */
void
copy_props( struct measurements *M , 
	    const size_t Nprops ) ;

/**
   @fn int compute_correlator( struct measurements *M , const size_t stride1 , const size_t stride2 , const size_t tshifted )
   @brief compute the momentum-projected correlation function in @M.corr
 */
int
compute_correlator( struct measurements *M , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted ) ;

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
