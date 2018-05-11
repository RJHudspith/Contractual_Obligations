#ifndef SPIN_INDEPENDENT_H
#define SPIN_INDEPENDENT_H

void
term_C0( struct NRQCD_fields *F ,
	 const size_t t ,
	 struct NRQCD_params NRQCD ) ;

void
term_C1_C6( struct NRQCD_fields *F ,
	    const size_t t ,
	    struct NRQCD_params NRQCD ) ;

void
term_C2( struct NRQCD_fields *F ,
	 const size_t t ,
	 struct NRQCD_params NRQCD ) ;

void
term_C5( struct NRQCD_fields *F ,
	 const size_t t ,
	 struct NRQCD_params NRQCD ) ;

void
term_C10EB( struct NRQCD_fields *F ,
	    const size_t t ,
	    struct NRQCD_params NRQCD ) ;

void
term_C11( struct NRQCD_fields *F ,
	  const size_t t ,
	  struct NRQCD_params NRQCD ) ;

#endif
