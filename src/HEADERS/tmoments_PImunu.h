/**
   @file tmoments_PImunu.h
   @brief prototype functions for temporal correlators
 */

#ifndef TMOMENTS_PIMUNU_H
#define TMOMENTS_PIMUNU_H

/**
   @fn void tmoments( const struct PIdata *data , const projtype PROJ ) 
   @brief time moments computation
 */
void
tmoments( const struct PIdata *data ,
	  const projtype PROJ ) ;

#endif