#ifndef GRAD_2_H
#define GRAD_2_H

void
gradsq( struct halfspinor *der2 ,
	const struct halfspinor *S ,
	const size_t i ,
	const size_t t ) ;

void
grad_sq_LCU( struct halfspinor *der2 ,
	     const struct halfspinor *S ,
	     const size_t t ) ;

void
gradsq_imp( struct halfspinor *der ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ) ;

#endif
