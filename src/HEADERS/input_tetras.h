/**
   @file input_tetras.h
   @brief prototype declarations for tetraquark contraction information
 */
#ifndef INPUT_TETRAS_H
#define INPUT_TETRAS_H

/**
   @fn int tetra_contractions( struct tetra_info *tetras , int *ntetras , struct inputs *INPUT , const int nprops , const GLU_bool first_pass )
   @brief get the tetraquark information from the input file
   @return #SUCCESS or #FAILURE
 */
int
tetra_contractions( struct tetra_info *tetras , 
		    int *ntetras ,
		    struct inputs *INPUT ,
		    const int nprops ,
		    const GLU_bool first_pass ) ;

#endif
