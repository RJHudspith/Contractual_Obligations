/**
   @file definitions.h
   @brief some useful macro definitions
 */
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/**
   @def B_CHANNELS
   @brief number of baryon channels we look at
 */
#define B_CHANNELS 6

/**
   @def CONDOR_MODE
   @brief if we are distributing this over multiple architectures
 */
#ifndef NOT_CONDOR_MODE
  #define CONDOR_MODE
#endif

/**
   @def FAILURE
   @brief our flag for when shit hits the fan
 */
#define FAILURE -1

/**
   @def GLU_STR_LENGTH
   @brief max length of a token in the input file
 */
#define GLU_STR_LENGTH 256

/**
   @def LT
   @brief Length of the time direction
 */
#define LT ( Latt.dims[ ND - 1 ] )

/**
   @def LCU
   @brief spatial volume
 */
#define LCU VOL3

/**
   @def MAX_CONTRACTIONS
   @brief break for the while loops in contraction logic
*/
#define MAX_CONTRACTIONS 20

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
   @def LVOLUME
   @brief lattice volume
 */
#define LVOLUME VOL4

/**
   @def NC
   @brief number of colors in our theory
 */
#ifndef NC
  #define NC 3
#endif

/**
   @def NCNC
   @brief Jamie stores links in a flattened #NC*#NC array
 */
#define NCNC ( NC*NC )

/**
   @def ND
   @brief number of dimensions of our theory
 */
#ifndef ND
  #define ND 4
#endif

/**
   @def NS
   @brief number of spins in our theory
 */
#ifndef NS
  #define NS 4
#endif

/**
   @def NSNS
   @brief sometimes we want to look at the whole NS*NS component of a spinor
 */
#define NSNS ( NS*NS )

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

/**
   @def SUCCESS
   @brief anything that isn't a failure is a success in our eyes
 */
#define SUCCESS !FAILURE

/**
   @def TETRA_NOPS
   @brief number of tetraquark operators we consider
 */
#define TETRA_NOPS 5

/**
   @def TWOPI
   @brief \f$ 2 \pi$ appears everywhere
 */
#define TWOPI 6.283185307179586

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
   @def VPF_MAGIC
   @brief magic number for this library's data
   spells VPF in ASCII
 */
#define VPF_MAGIC 717685

#endif
