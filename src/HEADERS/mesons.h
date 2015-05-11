/**
   @file mesons.h
   @brief prototype declerations for dispersion relation calculation
 */
#ifndef MESONS_H
#define MESONS_H

/**
   @fn int mesons_diagonal( struct propagator prop , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour diagonal meson dispersion relation
 */
int
mesons_diagonal( struct propagator prop ,
		 const struct cut_info CUTINFO ,
		 const char *outfile ) ;

#endif
