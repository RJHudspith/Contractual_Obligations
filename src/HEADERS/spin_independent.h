#ifndef SPIN_INDEPENDENT_H
#define SPIN_INDEPENDENT_H

void
term_C0( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C1_C6( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD ) ;

void
term_C2( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C5( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C10EB( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD ) ;

void
term_C11( struct NRQCD_fields *F ,
	  const size_t t ,
	  struct NRQCD_params NRQCD ) ;

#endif
