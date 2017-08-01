/**
   @file penta_udusb.h
   @brief prototype function declarations for a udusb - type pentaquark
 */
#ifndef PENTA_UDUSB_H
#define PENTA_UDUSB_H

/**
   @fn int pentaquark_udusb( struct propagator prop1 , struct propagator prop2 , struct propagator prop3 , struct cut_info CUTINFO , const char *outfile )
   @brief pentaquark candidate contraction
   @param prop1 :: light quark propagator
   @param prop2 :: strange quark propagator
   @param prop3 :: heavy quark propagator
   @param CUTINFO :: information on momentum cuts
   @param outfile :: output file name

   Contracts a ( u d )( u s ) \bar{b} - type pentaquark
 */
int
pentaquark_udusb( struct propagator prop1 , // L
		  struct propagator prop2 , // S
		  struct propagator prop3 , // H
		  struct cut_info CUTINFO ,
		  const char *outfile ) ;

#endif
