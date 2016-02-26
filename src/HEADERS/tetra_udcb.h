/**
   @file tetra_udcb.h
   @brief tetraquark contraction where two light quarks are degenerate but the heavy quarks aren't
 */
#ifndef TETRA_UDCB_H
#define TETRA_UDCB_H

/**
   @fn int tetraquark_udcb( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction where light content is degenerate
   @return #SUCCESS or #FAILURE
 */
int
tetraquark_udcb( struct propagator prop1 , // L1
		 struct propagator prop2 , // H1
		 struct propagator prop3 , // H2
		 struct cut_info CUTINFO ,
		 const char *outfile ) ;

#endif
