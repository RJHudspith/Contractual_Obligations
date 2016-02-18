/**
   @file diquark.h
   @brief diquark-diquark contraction
 */
#ifndef DIQUARK_H
#define DIQUARK_H

/**
   @fn int diquark( struct propagator S1 , struct propagator S2 , const struct cut_info CUTINFO , const char *outfile )
   @brief non-degenerate diquark contraction
 */
int
diquark( struct propagator S1 ,
	 struct propagator S2 ,
	 const struct cut_info CUTINFO , 
	 const char *outfile ) ;

#endif
