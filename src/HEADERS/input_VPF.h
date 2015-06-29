/**
   @file input_VPF.h
   @brief prototype declarations for VPF contraction logic
 */
#ifndef INPUT_VPF_H
#define INPUT_VPF_H

/**
   @fn int VPF_contractions( struct VPF_info *VPF , int *nVPF , const struct inputs *INPUT , const int nprops , const GLU_bool first_pass ) 
   @brief VPF contraction logic
 */
int
VPF_contractions( struct VPF_info *VPF , 
		  int *nVPF ,
		  const struct inputs *INPUT ,
		  const int nprops ,
		  const GLU_bool first_pass ) ;

#endif
