/**
   @file bar_ops_SSE.h
   @brief prototype function declarations for SSEd baryon operations
 */
#ifndef BAR_OPS_SSE_H
#define BAR_OPS_SSE_H

/**
   @fn const double complex baryon_contract( const struct spinor DiQ , const struct spinor S , const int d0 , const int d1 , const int d2 , const int d3 )
   @brief color trace a diquark with a propagator
   @return the color trace
 */
// otherwise do the usual
const double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const int d0 ,
		 const int d1 ,
		 const int d2 ,
		 const int d3 ) ;

/**
   @fn void cross_color_trace( struct spinor *__restrict DiQ , const struct spinor S ) 
   @brief color cross product and writes back into the Diquark
 */
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S ) ;

#endif
