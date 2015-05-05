/**
   @file input_disprel.h
   @brief prototype declarations for dispersion relation measurements
 */
#ifndef INPUT_DISPREL_H
#define INPUT_DISPREL_H

/**
   @fn int dispersion_contractions( struct dispersion_info *dispersions , int *ndispersions , const struct inputs *INPUT , const int nprops , const GLU_bool first_pass ) 
   @brief dispersion relation contractions
   @return #SUCCESS or #FAILURE
 */
int
dispersion_contractions( struct dispersion_info *dispersions , 
			 int *ndispersions ,
			 const struct inputs *INPUT ,
			 const int nprops ,
			 const GLU_bool first_pass ) ;

#endif
