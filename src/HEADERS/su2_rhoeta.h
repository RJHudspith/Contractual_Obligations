/**
   @file su2_rhoeta.h
   @brief prototype declarations for su2_rhoeta calculation
 */
#ifndef SU2_RHOETA_H
#define SU2_RHOETA_H

/**
   @fn int su2_rhoeta( struct propagator prop1 , struct cut_info CUTINFO , const char *outfile )
   @brief perform the rho-eta su2 contractions
 */
int
su2_rhoeta( struct propagator prop1 ,
	    struct cut_info CUTINFO ,
	    const char *outfile ) ;

#endif
