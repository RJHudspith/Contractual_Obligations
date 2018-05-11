#ifndef DERIVS_H
#define DERIVS_H

void
grad_sq( struct halfspinor *der2 ,
	 const struct halfspinor *S ,
	 const size_t t ) ;
        
void
grad2( struct halfspinor *der2 ,
       const struct halfspinor *S ,
       const size_t t ,
       const size_t mu ) ;

void
grad( struct halfspinor *der ,
      const struct halfspinor *S ,
      const size_t t ,
      const size_t mu ) ;

void
grad_imp( struct halfspinor *der ,
	  const struct halfspinor *S ,
	  const size_t t ,
	  const size_t mu ) ;

void
grad_sq_imp( struct halfspinor *der ,
	     const struct halfspinor *S ,
	     const size_t t ) ;

#endif
