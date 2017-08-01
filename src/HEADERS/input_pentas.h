/**
   @file input_pentas.h
   @brief prototype declarations for reading penta inputs
 */
#ifndef INPUT_PENTAS_H
#define INPUT_PENTAS_H

/**
   @fn int penta_contractions( struct tetra_info *pentas , size_t *npentas , struct inputs *INPUT , const size_t nprops , const GLU_bool first_pass ) 
   @brief read the input file to get the pentaquark contractions
   @param pentas :: pentaquark contraction map
   @param npentas :: number of pentaquark contractions
   @param INPUT :: input file data
   @param nprops :: number of propagators listed in the input file
   @param first_pass :: is this the first pass of reading the input file?
 */
int
penta_contractions( struct penta_info *pentas , 
		    size_t *npentas ,
		    struct inputs *INPUT ,
		    const size_t nprops ,
		    const GLU_bool first_pass ) ;

#endif
