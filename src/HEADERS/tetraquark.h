/**
   @file tetraquark.h
   @brief prototype tetraquark functions
 */
#ifndef TETRAQUARK_H
#define TETRAQUARK_H

/**
   @fn int tetraquark( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , const struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction code
   @return #SUCCESS or #FAILURE
 */
int
tetraquark( struct propagator prop1 ,
	    struct propagator prop2 ,
	    struct propagator prop3 ,
	    const struct cut_info CUTINFO ,
	    const char *outfile ) ;

#endif
