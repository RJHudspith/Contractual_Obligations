/**
   @file penta_bubds.h
   @brief prototype function declarations for a bubds - type pentaquark
 */
#ifndef PENTA_BUBDS_H
#define PENTA_BUBDS_H

/**
   @fn int pentaquark_bubds( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , struct cut_info CUTINFO , const char *outfile )
   @brief pentaquark candidate contraction
   @param prop1 :: heavy quark propagator
   @param prop2 :: light quark propagator
   @param prop3 :: some other anti quark propagator
   @param CUTINFO :: information on momentum cuts
   @param outfile :: output file name

   Contracts a ( b u )( b d ) \bar{s} - type pentaquark
 */
int
pentaquark_bubds( struct propagator prop1 , // H
		  struct propagator prop2 , // L
		  struct propagator prop3 , // whatever
		  struct cut_info CUTINFO ,
		  const char *outfile ) ;

#endif
