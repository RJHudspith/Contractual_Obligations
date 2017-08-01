/**
   @file wrap_pentas.h
   @brief pentaquark contraction logic prototype declarations
 */
#ifndef WRAP_PENTAS_H
#define WRAP_PENTAS_H

/**
   @fn int contract_pentas( struct propagator *prop , const struct penta_info *pentas , const struct cut_info CUTINFO , const size_t npentas )
   @brief contraction of pentaquark candidates
   @param prop :: propagator array
   @param pentas :: pentaquark contraction information
   @param CUTINFO :: momentum space cutting information
   @param npentas :: number of pentaquark contractions to do
 */
int
contract_pentas( struct propagator *prop ,
		 const struct penta_info *pentas ,
		 const struct cut_info CUTINFO ,
		 const size_t npentas ) ;

#endif
