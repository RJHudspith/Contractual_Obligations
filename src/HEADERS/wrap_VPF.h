/**
   @file wrap_VPF.h
   @brief wrapper for the various VPF contraction types
 */
#ifndef WRAP_VPF_H
#define WRAP_VPF_H

/**
   @fn int contract_VPF( struct propagator *prop , const struct site *lat , const struct VPF_info *VPF , const size_t nVPF , const struct cut_info CUTINFO )
   @brief wrapper for performing VPF contractions
   @return #SUCCESS or #FAILURE
 */
int
contract_VPF( struct propagator *prop ,
	      const struct site *lat ,
	      const struct VPF_info *VPF ,
	      const struct cut_info CUTINFO ,
	      const size_t nVPF ) ;

#endif
