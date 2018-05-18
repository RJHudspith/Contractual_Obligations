#ifndef GRAD_4_H
#define GRAD_4_H

void
grad4( struct halfspinor *der ,
       const struct halfspinor *S ,
       const struct field *Fmunu ,
       const size_t i ,
       const size_t t ,
       const size_t mu ) ;

void
grad_sqsq( struct halfspinor *der2 ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t ) ;

#endif
