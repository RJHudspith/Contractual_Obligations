/**
   @file definitions.h
   @brief some useful macro definitions
 */
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/**
   @def FAILURE
   @brief our flag for when shit hits the fan
 */
#define FAILURE -1

/**
   @def SUCCESS
   @brief anything that isn't a failure is a success in our eyes
 */
#define SUCCESS !FAILURE

/**
   @def MAX_LINE_LENGTH
   @brief maximun NERSC header token length
 */
#define MAX_LINE_LENGTH 64

/**
   @def MAX_TOKENS
   @brief maximum number of NERSC header values
 */
#define MAX_TOKENS 32

/**
   @def GLU_STR_LENGTH
   @brief max length of a token in the input file
 */
#define GLU_STR_LENGTH 128

/**
   @def L0
   @brief Length of the time direction
 */
#define L0 ( Latt.dims[ ND - 1 ] )

/**
   @def VOL3 
   @brief spatial hypercube size
 */
#define VOL3 Latt.Lcu

/**
   @def VOL4
   @brief full lattice volume
 */
#define VOL4 Latt.Volume

/**
   @def LCU
   @brief spatial volume
 */
#define LCU VOL3

/**
   @def LVOLUME
   @brief lattice volume
 */
#define LVOLUME VOL4

/**
   @def TWOPI
   @brief \f$ 2 \pi$ appears everywhere
 */
#define TWOPI ( 2.0 * M_PI )

#ifndef NC
/**
   @def NC
   @brief number of colors in our theory
 */
  #define NC 3
#endif

#ifndef NS
/**
   @def NS
   @brief number of spins in our theory
 */
  #define NS 4
#endif

/**
   @def ND
   @brief number of dimensions of our theory
 */
#ifndef ND
  #define ND 4
#endif

/**
   @def NSNS
   @brief sometimes we want to look at the whole NS*NS component of a spinor
 */
#define NSNS ( NS*NS )

/**
   @def NCNC
   @brief Jamie stores links in a flattened #NC*#NC array
 */
#define NCNC ( NC*NC )

/**
   @def PLAQ_AND_TRACE_TOL
   @brief used to check the correctness of the NERSC configuration
 */
#define PLAQ_AND_TRACE_TOL 1E-6

/**
   @def PREC_TOL
   @brief precision tolerance 
 */
#define PREC_TOL (NC * 1.0E-14)

#endif
