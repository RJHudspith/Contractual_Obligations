/**
   @file input_VPF.h
   @brief prototype declarations for VPF contraction logic
 */
#ifndef INPUT_VPF_H
#define INPUT_VPF_H

/**
   @fn int VPF_contractions( struct VPF_info *VPF , size_t *nVPF , const struct inputs *INPUT , const size_t nprops , const GLU_bool first_pass ) 
   @brief VPF contraction logic
 */
int
VPF_contractions( struct VPF_info *VPF , 
		  size_t *nVPF ,
		  const struct inputs *INPUT ,
		  const size_t nprops ,
		  const GLU_bool first_pass ) ;

#endif
