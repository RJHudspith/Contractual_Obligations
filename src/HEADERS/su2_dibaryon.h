/**
   @file su2_dibaryon.h
   @brief prototype declarations for su2 dibaryon contraction
 */
#ifndef SU2_DIBARYON_H
#define SU2_DIBARYON_H

/**
   @fn int su2_dibaryon( struct propagator prop1 , struct cut_info CUTINFO , const char *outfile )
   @brief su2 dibaryon contraction code
 */
int
su2_dibaryon( struct propagator prop1 ,
	      struct cut_info CUTINFO ,
	      const char *outfile ) ;

#endif
