/**
   @file tetra_uscb.h
   @brief tetraquark contraction where all quarks are non-degenerate
 */
#ifndef TETRA_USCB_H
#define TETRA_USCB_H

/**
   @fn int tetraquark_uscb( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , struct propagator prop4 , struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction where light content is degenerate
   @return #SUCCESS or #FAILURE
 */
int
tetraquark_uscb( struct propagator prop1 ,
		 struct propagator prop2 ,
		 struct propagator prop3 ,
		 struct propagator prop4 ,
		 struct cut_info CUTINFO ,
		 const char *outfile ) ;

#endif
