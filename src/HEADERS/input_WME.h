/**
   @file input_WME.h
   @brief prototype declarations for WME contraction logic
 */
#ifndef INPUT_WME_H
#define INPUT_WME_H

/**
   @fn int matrix_element_contractions( struct WME_info *wme , int *nWME , const struct inputs *INPUT , const int nprops , const GLU_bool first_pass ) 
   @brief WME contractions
   @return #SUCCESS or #FAILURE
 */
int
matrix_element_contractions( struct WME_info *wme , 
			     int *nWME ,
			     const struct inputs *INPUT ,
			     const int nprops ,
			     const GLU_bool first_pass ) ;

#endif
