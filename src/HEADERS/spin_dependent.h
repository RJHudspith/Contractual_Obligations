#ifndef SPIN_DEPENDENT_H
#define SPIN_DEPENDENT_H

void
term_C3( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

void
term_C4( struct NRQCD_fields *F ,
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

#endif
