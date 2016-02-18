/**
   @file diquark_degen.h
   @brief diquark-diquark contraction prototypes
 */
#ifndef DIQUARK_DEGEN_H
#define DIQUARK_DEGEN_H

/**
   @fn int diquark_degen( struct propagator prop1 , const struct cut_info CUTINFO , const char *outfile )
   @brief diquark-diquark contraction where the quarks are degenerate
 */
int
diquark_degen( struct propagator S1 ,
	       const struct cut_info CUTINFO , 
	       const char *outfile ) ;

#endif
