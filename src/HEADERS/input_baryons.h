/**
   @file input_baryons.h
   @brief prototype declarations for baryon contraction instructions from the input file
 */
#ifndef INPUT_BARYONS_H
#define INPUT_BARYONS_H

/**
   @fn int baryon_contractions( struct baryon_info *baryons , size_t *nbaryons , const struct inputs *INPUT , const size_t nprops , const GLU_bool first_pass ) 
   @brief baryon contraction information
   @return #SUCCESS or #FAILURE
 */
int
baryon_contractions( struct baryon_info *baryons , 
		     size_t *nbaryons ,
		     const struct inputs *INPUT ,
		     const size_t nprops ,
		     const GLU_bool first_pass ) ;

#endif
