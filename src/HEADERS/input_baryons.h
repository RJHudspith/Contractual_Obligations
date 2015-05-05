/**
   @file input_baryons.h
   @brief prototype declarations for baryon contraction instructions from the input file
 */
#ifndef INPUT_BARYONS_H
#define INPUT_BARYONS_H

/**
   @fn int baryon_contractions( struct baryon_info *baryons , int *nbaryons , const struct inputs *INPUT , const int nprops , const GLU_bool first_pass ) 
   @brief baryon contraction information
   @return #SUCCESS or #FAILURE
 */
int
baryon_contractions( struct baryon_info *baryons , 
		     int *nbaryons ,
		     const struct inputs *INPUT ,
		     const int nprops ,
		     const GLU_bool first_pass ) ;

#endif
