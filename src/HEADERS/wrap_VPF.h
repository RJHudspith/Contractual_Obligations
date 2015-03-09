/**
   @file wrap_VPF.h
   @brief wrapper for the various VPF contraction types
 */

#ifndef WRAP_VPF_H
#define WRAP_VPF_H

/**
   @fn int contract_VPF( FILE **fprops , const struct site *lat , const struct VPF_info *VPF , const int nVPF , const struct cut_info CUTINFO )
   @brief wrapper for performing VPF contractions

int contract_VPF( FILE **fprops , const struct meson_info *mesons , const int nmesons )
   @return #SUCCESS or #FAILURE
 */
int
contract_VPF( FILE **fprops ,
	      const struct site *lat ,
	      const struct VPF_info *VPF ,
	      const int nVPF ,
	      const struct cut_info CUTINFO ) ;

#endif
