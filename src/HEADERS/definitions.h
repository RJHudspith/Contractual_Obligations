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
#ifndef B_CHANNELS
  #define B_CHANNELS (16)
#endif
// if it is defined to something silly we guard against that also
#if (B_CHANNELS < 0) || (B_CHANNELS>16)
  #define B_CHANNELS (16)
#endif

/**
   @def M_CHANNELS
   @brief number of meson channels we look at
 */
#ifndef M_CHANNELS
  #define M_CHANNELS (16)
#endif
// if it is defined to something silly we guard against that also
#if (M_CHANNELS < 0) || (M_CHANNELS>16)
  #define M_CHANNELS (16)
#endif

/**
   @def MWC_4096_RNG
   @brief the type of parallel RNG we use
 */
#define MWC_4096_RNG

/**
   @def RNG_TABLE
   @brief length of the rng table available
 */
#define RNG_TABLE 4096

/**
   @def CONDOR_MODE
   @brief if we are distributing this over multiple architectures
 */
#ifndef NOT_CONDOR_MODE
  #define CONDOR_MODE
#endif

/**
   @def CORR_MAGIC
   @brief CORR magic number for correlator files
 */
#define CORR_MAGIC (67678282)

/**
   @def FAILURE
   @brief our flag for when shit hits the fan
 */
#define FAILURE -1

/**
   @def GLU_STR_LENGTH
   @brief max length of a token in the input file
 */
#define GLU_STR_LENGTH (256)

/**
   @def IO_NBLOCK
   @brief size of blocking factor in IO
 */
#ifndef IO_NBLOCK
  #define IO_NBLOCK (1)
#endif

/**
   @def LT
   @brief Length of the time direction
 */
#define LT (Latt.dims[ND-1])

/**
   @def LCU
   @brief spatial volume
 */
#define LCU (Latt.Lcu)

/**
   @def LVOLUME
   @brief lattice volume
 */
#define LVOLUME (Latt.Volume)

/**
   @def MAX_CONTRACTIONS
   @brief break for the while loops in contraction logic
*/
#define MAX_CONTRACTIONS (64)

/**
   @def MAX_LINE_LENGTH
   @brief maximun NERSC header token length
 */
#define MAX_LINE_LENGTH (64)

/**
   @def MAX_TOKENS
   @brief maximum number of NERSC header values
 */
#define MAX_TOKENS (64)

/**
   @def NC
   @brief number of colors in our theory
 */
#ifndef NC
  #define NC (3)
#endif

/**
   @def NCNC
   @brief Jamie stores links in a flattened #NC*#NC array
 */
#define NCNC (NC*NC)

/**
   @def ND
   @brief number of dimensions of our theory
 */
#ifndef ND
  #define ND (4)
#endif

/**
   @def NRQCD_TOL
   @brief tolerance at which we consider NRQCD parameters to be zero
 */
#define NRQCD_TOL (1.E-6)

/**
   @def NS
   @brief number of spins in our theory
 */
#ifndef NS
  #define NS (4)
#endif

/**
   @def NSNS
   @brief sometimes we want to look at the whole NS*NS component of a spinor
 */
#define NSNS (NS*NS)

/**
   @def PLAQ_AND_TRACE_TOL
   @brief used to check the correctness of the NERSC configuration
 */
#define PLAQ_AND_TRACE_TOL (1E-6)

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
   @def PENTA_NOPS
   @brief number of pentaquark operators
 */
#define PENTA_NOPS (9)

/**
   @def PENTA_NBLOCK
   @brief number of different gamma combinations for the penta operator
 */
#define PENTA_NBLOCK (2)

/**
   @def TETRA_NBLOCK
   @brief number of different tetraquark blocks we will use
 */
#ifndef TETRA_NBLOCK
  #define TETRA_NBLOCK (1)
#endif

/**
   @def TETRA_NOPS
   @brief number of tetraquark operators (TETRA_NBLOCK)^2x4 we consider
 */
#define TETRA_NOPS (TETRA_NBLOCK*TETRA_NBLOCK*4)

/**
   @def TWOPI
   @brief \f$ 2 \pi$ appears everywhere
 */
#define TWOPI (6.283185307179586)

/**
   @def VPF_MAGIC
   @brief magic number for this library's data
   spells VPF in ASCII
 */
#define VPF_MAGIC (717685)

/**
   @def PENTA_NCOLORS
   @brief number of colour indices we use in the penta contractions
 **/
#define PENTA_NCOLORS (6561)
  
#endif
