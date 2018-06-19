/**
   @file sources.h
   @brief prototype declarations for implementening NRQCD sources
*/
#ifndef SOURCES_H
#define SOURCES_H

/**
   @fn void initialise_source( struct halfspinor *S , const struct propagator prop )
   @brief initialise a source into S
 */
void
initialise_source( struct halfspinor *S ,
		   const struct propagator prop ) ;

#endif
