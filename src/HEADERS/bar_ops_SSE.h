/**
   @file bar_ops_SSE.h
   @brief prototype function declarations for SSEd baryon operations
 */
#ifndef BAR_OPS_SSE_H
#define BAR_OPS_SSE_H

/**
   @fn double complex baryon_contract( const struct spinor DiQ , const struct spinor S , const size_t d0 , const size_t d1 , const size_t d2 , const size_t d3 )
   @brief color trace a diquark with a propagator
   @return the color trace
 */
// otherwise do the usual
double complex
baryon_contract( const struct spinor DiQ ,
		 const struct spinor S ,
		 const size_t d0 ,
		 const size_t d1 ,
		 const size_t d2 ,
		 const size_t d3 ) ;

/**
   @fn void cross_color_trace( struct spinor *__restrict DiQ , const struct spinor S ) 
   @brief color cross product and writes back into the Diquark
 */
void
cross_color_trace( struct spinor *__restrict DiQ ,
		   const struct spinor S ) ;

#endif
