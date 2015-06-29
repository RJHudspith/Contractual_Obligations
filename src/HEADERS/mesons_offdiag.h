/**
   @file mesons_offdiag.h
   @brief prototype declarations for flavour off-diagonal mesons
 */
#ifndef MESONS_OFFDIAG_H
#define MESONS_OFFDIAG_H

/**
   @fn int mesons_offdiagonal( struct propagator prop1 , struct propagator prop2 , const struct cut_info CUTINFO , const char *outfile )
   @brief flavour off-diagonal meson dispersion relation
 */
int
mesons_offdiagonal( struct propagator prop1 ,
		    struct propagator prop2 ,
		    const struct cut_info CUTINFO ,
		    const char *outfile ) ;

#endif
