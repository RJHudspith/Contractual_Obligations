#ifndef DERIVS_H
#define DERIVS_H

void
grad_imp_LCU( struct halfspinor *der ,
	      const struct halfspinor *S ,
	      const struct field *Fmunu ,
	      const size_t t ,
	      const size_t mu ) ;

void
FMUNU_grad_imp( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx ) ;

void
grad_imp_FMUNU( struct halfspinor *der ,
		const struct halfspinor *S ,
		const struct field *Fmunu ,
		const size_t i ,
		const size_t t ,
		const size_t mu ,
		const size_t Fmunu_idx ) ;

#endif
