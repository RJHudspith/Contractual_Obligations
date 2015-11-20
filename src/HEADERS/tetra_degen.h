/**
   @file tetra_degen.h
   @brief tetraquark contraction where two light quarks are degenerate
 */
#ifndef TETRA_DEGEN_H
#define TETRA_DEGEN_H

/**
   @fn int tetraquark_degen( struct propagator prop1 , struct propagator prop2 , struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction where light content is degenerate
   @return #SUCCESS or #FAILURE
 */
int
tetraquark_degen( struct propagator prop1 ,
		  struct propagator prop2 ,
		  struct cut_info CUTINFO ,
		  const char *outfile ) ;

#endif
