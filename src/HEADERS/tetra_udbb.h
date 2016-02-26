/**
   @file tetra_udbb.h
   @brief tetraquark contraction where two light quarks are degenerate as are the two heavy quarks
 */
#ifndef TETRA_UDBB_H
#define TETRA_UDBB_H

/**
   @fn int tetraquark_udbb( struct propagator prop1 , struct propagator prop2 , struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction where light content is degenerate
   @return #SUCCESS or #FAILURE
 */
int
tetraquark_udbb( struct propagator prop1 ,
		  struct propagator prop2 ,
		  struct cut_info CUTINFO ,
		  const char *outfile ) ;

#endif
