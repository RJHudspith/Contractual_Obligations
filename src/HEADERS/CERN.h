/**
   @file CERN.h
   @brief CERN configuration file reader
 */
#ifndef CERN_H
#define CERN_H

/**
   @fn int read_CLS_field( struct site *__restrict lat , FILE *__restrict in , uint32_t *chksum )
   @brief read a CLS/CERN configuration file
 */
int
read_CLS_field( struct site *__restrict lat , 
		FILE *__restrict in , 
		uint32_t *chksum ) ;

#endif
