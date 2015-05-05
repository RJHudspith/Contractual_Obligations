/**
   @file dispersions.h
   @brief prototype declerations for dispersion relation calculation
 */
#ifndef DISPERSIONS_H
#define DISPERSIONS_H

/**
   @fn int dispersions_diagonal( struct propagator prop , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour diagonal meson dispersion relation
 */
int
dispersions_diagonal( struct propagator prop ,
		      const struct cut_info CUTINFO ,
		      const char *outfile ) ;

/**
   @fn int dispersions_offdiagonal( struct propagator prop1 , struct propagator prop2 , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour off-diagonal meson dispersion relation
 */
int
dispersions_offdiagonal( struct propagator prop1 ,
			 struct propagator prop2 ,
			 const struct cut_info CUTINFO ,
			 const char *outfile ) ;

#endif
