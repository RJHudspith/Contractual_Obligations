/**
   @file wrap_tetras.h
   @brief prototype functions for tetraquark logic
 */
#ifndef WRAP_TETRAS_H
#define WRAP_TETRAS_H

/**
   @fn int contract_tetras( struct propagator *prop , const struct tetra_info *tetras , const struct cut_info CUTINFO , const int ntetras )
   @brief tetraquark logic wrapper
 */
int
contract_tetras( struct propagator *prop ,
		 const struct tetra_info *tetras ,
		 const struct cut_info CUTINFO ,
		 const int ntetras ) ;

#endif
