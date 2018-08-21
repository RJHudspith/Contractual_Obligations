/**
   @file sources.h
   @brief prototype declarations for implementening NRQCD sources
*/
#ifndef SOURCES_H
#define SOURCES_H

/**
   @fn int initialise_source( struct halfspinor *S , struct halfspinor *S1 , const struct propagator prop )
   @brief initialise a source into S
   @return #SUCCESS or #FAILURE
 */
int
initialise_source( struct halfspinor *S ,
		   struct halfspinor *S1 ,
		   const struct propagator prop ) ;

#endif
