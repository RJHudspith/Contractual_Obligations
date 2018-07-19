/**
   @file spin_independent.h
   @brief prototype declarations for spin-dependent NRQCD terms
 */
#ifndef SPIN_INDEPENDENT_H
#define SPIN_INDEPENDENT_H

/**
   @fn void term_C0( struct halfspinor *H , const struct halfspinor *S , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_0 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"
   Computes \f$ -\frac{c_0}{2M_0}\Delta^{(2)} $\f
 */
void
term_C0( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C1_C6( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_1 & c_6 terms of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ -\left( \frac{c_1}{(2M_0)^3} + \frac{c_6}{4n(2M_0)^2} \right) (\tilde\Delta^{(2)})^2 $\f
 */
void
term_C1_C6( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C2( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_2 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ +i\frac{c_2}{2(2M_0)^2}  (\tilde\Delta\cdot\tilde{E} - \tilde{E} \cdot\tilde\Delta) $\f
 */
void
term_C2( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C5( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_5 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ +\frac{1}{24M_0}\tilde\Delta^{(4)} $\f
 */
void
term_C5( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C10EB( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_10EB term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"
   Computes \f$ -\frac{c_10}{(2M_0)^3}\left( \tilde{E}\cdot\tilde{E} + \tilde{B}\cdot\tilde{B} \right) $\f

   @warning does not split E.E and B.B terms. Needs vacuum subtraction.
 */
void
term_C10EB( struct halfspinor *H ,
	    const struct halfspinor *S ,
	    const struct field *Fmunu ,
	    const size_t i ,
	    const size_t t ,
	    const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C11( struct NRQCD_fields *F , const size_t t , struct NRQCD_params NRQCD )
   @brief c_11 term of NRQCD hamiltonian evaluated for all sites on timeslice "t"
   Computes \f$ -\frac{1}{24(n)^2(2M_0)^3}((\tilde\Delta^{(2)})^2)^2 $\f
 */
void
term_C11( struct NRQCD_fields *F ,
	  const size_t t ,
	  struct NRQCD_params NRQCD ) ;

#endif
