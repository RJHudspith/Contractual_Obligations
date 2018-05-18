#ifndef SPIN_DEPENDENT_H
#define SPIN_DEPENDENT_H

void
term_C3( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C4( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C7( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C8( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C9EB( struct halfspinor *H ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD ) ;

#endif
