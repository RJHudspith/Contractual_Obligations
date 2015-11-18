/**
   @file tetra_contractions.h
   @brief prototype declarations for tetraquark contractions
 */
#ifndef TETRA_CONTRACTIONS_H
#define TETRA_CONTRACTIONS_H

/**
   @fn double complex diquark_diquark( const struct spinor U , const struct spinor D , const struct spinor B , const struct gamma *GAMMAS , const size_t mu )
   @brief compute a diquark-diquark system
   @return the trace of the two portions
 */
double complex
diquark_diquark( const struct spinor U ,
		 const struct spinor D ,
		 const struct spinor B , // full adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu ) ;

/**
   @fn double complex dimeson_dimeson( const struct spinor U , const struct spinor D , const struct spinor B , const struct gamma *GAMMAS , const size_t mu )
   @brief computes a meson-meson system
   @return the trace of the various propagators
 */
double complex
dimeson_dimeson( const struct spinor U ,  // u prop
		 const struct spinor D ,  // d prop
		 const struct spinor B , // adjoint of B
		 const struct gamma *GAMMAS ,
		 const size_t mu ) ;

#endif
