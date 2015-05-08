/**
   @file tetraquark.h
   @brief prototype tetraquark functions
 */
#ifndef TETRAQUARK_H
#define TETRAQUARK_H

/**
   @fn int tetraquark_diagonal( struct propagator prop , const struct cut_info CUTINFO , const char *outfile )
   @brief tetraquark contraction code
 */
int
tetraquark_diagonal( struct propagator prop ,
		     const struct cut_info CUTINFO ,
		     const char *outfile ) ;

#endif
