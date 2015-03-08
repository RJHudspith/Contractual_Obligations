/**
   @file enum.h
   @brief storage for all the enums
 */

#ifndef ENUM_H
#define ENUM_H

#define DEBUG

//#define CHROMA_DIRAC_CONVENTION

/**
   @enum gamma_labels
   @brief different conventions have different labelling
*/
#ifdef CHROMA_DIRAC_CONVENTION
// degrand rossi on a chroma basis
typedef enum {
  IDENTITY = 0 ,
  GAMMA_0  = 1 ,
  GAMMA_1  = 2 ,
  GAMMA_2  = 4 ,
  GAMMA_3  = 8 ,
  GAMMA_5  = 15 ,
} gamma_labels ;
#else
// Mainz labelling
typedef enum {
  GAMMA_0  = 0 ,
  GAMMA_1  = 1 ,
  GAMMA_2  = 2 ,
  GAMMA_3  = 3 ,
  IDENTITY = 4 ,
  GAMMA_5  = 5 ,
} gamma_labels ;
#endif

/**
   @enum sourcetype
   @brief propagator source type
 */
typedef enum {
  WALL , 
  POINT } sourcetype ;

/**
   @enum proptype
   @brief propagator type definition
 */
typedef enum {
  CHIRAL ,
  CHIRAL_TO_NREL ,
  NREL } proptype ;

/**
   @enum GLU_endian
   @brief enumeration of the endianess
*/
typedef enum 
  { L_ENDIAN , 
    B_ENDIAN } GLU_endian ;

/**
   @enum GLU_output
   @brief enumerate the output types
 */
typedef enum
  { 
    NO_OUTPUT ,
    OUTPUT_SMALL ,  // logically the fewest number of params saved 
    OUTPUT_GAUGE ,  // NERSC's gauge , the top NC-1 rows and then complete
    OUTPUT_NCxNC ,  // the whole matrix, wasteful
    OUTPUT_MILC ,   // the whole matrix, wasteful
    OUTPUT_ILDG ,   // the whole matrix, wasteful
    OUTPUT_SCIDAC , // the whole matrix, wasteful
    OUTPUT_HIREP } GLU_output ; // the whole matrix in HIREP's order

/**
   @enum file_prec
   @brief defs for the readers
*/
typedef enum 
  { FLOAT_PREC , 
    DOUBLE_PREC } file_prec ;

/**
   @enum GLU_bool
   @brief standard boolean type deal
   In C FALSE = 0 , TRUE = {ANYTHING_ELSE}
 */
typedef enum 
  { GLU_FALSE , 
    GLU_TRUE } GLU_bool ;

/**
   @enum header_mode
   @brief enums for the available headers my code can read
 */
typedef enum 
  { UNSUPPORTED ,  
    NERSC_HEADER ,
    HIREP_HEADER ,
    MILC_HEADER ,
    SCIDAC_HEADER ,
    ILDG_SCIDAC_HEADER ,
    ILDG_BQCD_HEADER ,
    LIME_HEADER ,
    RANDOM_CONFIG ,
    UNIT_GAUGE ,
    INSTANTON } header_mode ;

// enumerated type for the projection
typedef enum {
  TRANS_PROJ ,
  LONG_PROJ ,
  ZERO_PLUS_ONE } PImunu_projtype ;

// enumerated type for the WI correction index
typedef enum {
  CORR_MU ,
  CORR_NU ,
  CORR_MUpNU ,
  UNCORR } correction_dir ;

/**
   @enum momentum_cut_def
   @brief enumerate the cutting types
 */
typedef enum
  { HYPERCUBIC_CUT , 
    PSQ_CUT ,
    CYLINDER_CUT ,
    CYLINDER_AND_CONICAL_CUT } momentum_cut_def ;

#endif
