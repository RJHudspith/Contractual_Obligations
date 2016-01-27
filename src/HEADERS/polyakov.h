/**
   @file polyakov.h
   @brief polyakov loop (static quark) prototype function declarations
 */
#ifndef POLYAKOV_H
#define POLYAKOV_H

/**
   @fn static_quark( struct spinor *S , const struct gamma gamma_t , const size_t t , const GLU_bool forward )
   @brief computes a static quark propagator
 */
void
static_quark( struct spinor *S ,
	      const struct gamma gamma_t ,
	      const size_t t ,
	      const GLU_bool forward ) ;

/**
   @fn double complex poly( const struct site *__restrict lat , int dir )
   @brief trace of a polyakov loop in direction dir
 */
double complex
poly( const struct site *__restrict lat , 
      size_t dir ) ;

#endif
