/**
   @file quark_smear.h
   @brief quark source and sink smearing prototype declarations
 */
#ifndef QUARK_SMEAR_H
#define QUARK_SMEAR_H

/**
   @fn void sink_smear( struct spinor **S , struct spinor **S1 , const size_t t , const struct cut_info CUTINFO , const size_t Np )
   @brief perform sink smearing
 */
void
sink_smear( struct spinor **S ,
	    struct spinor **S1 ,
	    const size_t t ,
	    const struct cut_info CUTINFO ,
	    const size_t Np ) ;

/**
   @fn void source_smear( struct halfspinor *S , struct halfspinor *S1 , const size_t t , const struct source_info Source )
   @brief perform NRQCD source smearing
 */
void
source_smear( struct halfspinor *S ,
	      struct halfspinor *S1 ,
	      const size_t t ,
	      const struct source_info Source ) ;

#endif
