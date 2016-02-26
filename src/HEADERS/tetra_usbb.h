/**
   @file tetra_usbb.h
   @brief prototype tetraquark non-degenerate light, degenerate heavy
 */
#ifndef TETRA_USBB_H
#define TETRA_USBB_H

/**
   @fn int tetra_usbb( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , const struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction code
   @return #SUCCESS or #FAILURE
 */
int
tetraquark_usbb( struct propagator prop1 ,
		 struct propagator prop2 ,
		 struct propagator prop3 ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ) ;

#endif
