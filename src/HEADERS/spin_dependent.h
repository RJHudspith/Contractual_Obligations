/**
   @file spin_dependent.h
   @brief prototype declarations for spin-dependent NRQCD terms
 */
#ifndef SPIN_DEPENDENT_H
#define SPIN_DEPENDENT_H

/**
   @fn void term_C3( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_3 term of NRQCD Hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ -\frac{c_3}{2(2M_0)^2}\sigma\cdot\left( \tilde\Delta\times\tilde{E} - \tilde{E}\times\tilde\Delta \right)$\f where the tilde's mean O(a^2) improvement
 */
void
term_C3( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C4( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_4 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ -\frac{c_4}{2M_0}\sigma\cdot \tilde{B} $\f where the tilde implies O(a^2) improvement
 */
void
term_C4( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C7( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_7 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ -\frac{c_7}{(2M_0)^3}\left\{ \tilde\Delta^{(2)} , \sigma\cdot \tilde{B} \right\} $\f where the tilde implies O(a^2) improvement
 */
void
term_C7( struct halfspinor *H ,
	 const struct halfspinor *S ,
	 const struct field *Fmunu ,
	 const size_t i ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C8( struct NRQCD_fields *F , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_8 term of NRQCD hamiltonian evaluated for sites on timeslice "t"

   Computes \f$ -\frac{c_8}{4(2M_0)^4}\left\{ \tilde\Delta^{(2)} , \sigma\cdot\left( \tilde\Delta\times\tilde{E} - \tilde{E}\times\tilde\Delta \right) \right\}$\f where the tilde implies O(a^2) improvement
 */
void
term_C8( struct NRQCD_fields *F ,
	 const size_t t ,
	 const struct NRQCD_params NRQCD ) ;

/**
   @fn void term_C9EB( struct halfspinor *H , const struct halfspinor *S , const struct field *Fmunu , const size_t i , const size_t t , const struct NRQCD_params NRQCD )
   @brief c_9 term of NRQCD hamiltonian evaluated at site "i" on timeslice "t"

   Computes \f$ -\frac{c_9}{(2M_0)^3}\sigma\cdot\left( \tilde{E}\times\tilde{E} + \tilde{B}\times\tilde{B} \right)$\f where the tilde implies O(a^2) improvement

   @warning this does not treat E\times E and B\times B cases separately
 */
void
term_C9EB( struct halfspinor *H ,
	   const struct halfspinor *S ,
	   const struct field *Fmunu ,
	   const size_t i ,
	   const size_t t ,
	   const struct NRQCD_params NRQCD ) ;

#endif
