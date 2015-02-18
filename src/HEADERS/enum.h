/**
   @file enum.h
   @brief storage for all the enums
 */

#ifndef ENUM_H
#define ENUM_H

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

#endif