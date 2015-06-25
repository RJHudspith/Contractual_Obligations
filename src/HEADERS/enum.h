/**
   @file enum.h
   @brief storage for all the enums
 */
#ifndef ENUM_H
#define ENUM_H

/**
   @enum config_size
   @brief for writing out files
   loops number of spaces for this format uses NC and NCNC
 */
typedef enum
  { LOOP_SMALL = NCNC - 1 ,
    LOOP_GAUGE = 2 * ( NC - 1 ) * NC ,
    LOOP_NCxNC = 2 * NCNC } config_size ;

/**
   @enum correction_dir
   @brief type for the WI correction index
*/
typedef enum {
  CORR_MU ,
  CORR_NU ,
  CORR_MUpNU ,
  UNCORR } correction_dir ;

/**
   @enum current_type
   @brief fermionic current type
 */
typedef enum {
  LOCAL_LOCAL , 
  CONSERVED_LOCAL } current_type ;

/**
   @enum endianness
   @brief specify which end is up
 */
typedef enum
  { LILENDIAN , 
    BIGENDIAN } endianness ;

/**
   @enum file_prec
   @brief defs for the readers
*/
typedef enum 
  { FLOAT_PREC , 
    DOUBLE_PREC } file_prec ;

/**
   @enum fp_precision
   @brief specify the floating point precision
 */
typedef enum
  { SINGLE ,
    DOUBLE } fp_precision ;

/**
   @enum gamma_labels
   @brief different conventions have different labelling
*/
// Mainz labelling
typedef enum {
  GAMMA_0  = 0 ,
  GAMMA_1  = 1 ,
  GAMMA_2  = 2 ,
  GAMMA_3  = 3 ,
  IDENTITY = 4 ,
  GAMMA_5  = 5 } gamma_labels ;

/**
   @enum GLU_bool
   @brief standard boolean type deal
   In C FALSE = 0 , TRUE = {ANYTHING_ELSE}
 */
typedef enum 
  { GLU_FALSE , 
    GLU_TRUE } GLU_bool ;

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
   @enum header_mode
   @brief enums for the available headers my code can read
 */
typedef enum 
  { UNSUPPORTED = FAILURE ,  
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

/**
   @enum momentum_cut_def
   @brief enumerate the cutting types
 */
typedef enum
  { HYPERCUBIC_CUT , 
    PSQ_CUT ,
    CYLINDER_CUT ,
    CYLINDER_AND_CONICAL_CUT } momentum_cut_def ;

/**
   @enum PImunu_projtype
   @brief VPF projection types we support
 */
typedef enum {
  TRANS_PROJ ,
  LONG_PROJ ,
  ZERO_PLUS_ONE } PImunu_projtype ;

/**
   @enum proptype
   @brief propagator type definition
 */
typedef enum {
  CHIRAL ,
  NREL ,
  STATIC } proptype ;

/**
   @enum sourcetype
   @brief propagator source types recognised by us
 */
typedef enum
  { POINT ,
    WALL } sourcetype ;

/**
   @enum vector_axial
   @brief enumerate the currents we have available
 */
typedef enum
  { VECTOR , 
    AXIAL } vector_axial ;

#endif
