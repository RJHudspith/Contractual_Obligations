/**
   @file input_mesons.h
   @brief prototype declarations for meson contractions
 */
#ifndef INPUT_MESONS_H
#define INPUT_MESONS_H

/**
   @fn int meson_contractions( struct meson_info *mesons , int *nmesons , const struct inputs *INPUT , const int nprops , const GLU_bool first_pass )
   @brief meson contraction information, packs mesons struct
   @return #SUCCESS or #FAILURE
 */
int
meson_contractions( struct meson_info *mesons , 
		    int *nmesons ,
		    const struct inputs *INPUT ,
		    const int nprops ,
		    const GLU_bool first_pass ) ;

#endif
