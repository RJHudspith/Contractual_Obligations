#ifndef NRQCD_H
#define NRQCD_H

int
compute_nrqcd_props( struct propagator *prop ,
		     const size_t nprops ) ;

int
free_nrqcd_props( struct propagator *prop ,
		  const size_t nprops ) ;

#endif
