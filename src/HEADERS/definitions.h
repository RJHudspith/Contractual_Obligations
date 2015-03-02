/**
   @file definitions.h
   @brief some useful macro definitions
 */
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define FAILURE -1
#define SUCCESS !FAILURE

#define MAX_LINE_LENGTH 64

#define MAX_TOKENS 32

// we might want to change this at some point
#define GLU_STR_LENGTH 96

#define L0 ( Latt.dims[ ND - 1 ] )

#define VOL3 Latt.Lcu
#define VOL4 Latt.Volume

#define LCU Latt.Lcu
#define LVOLUME VOL4

#define TWOPI ( 2.0 * M_PI )

// default number of colors
#ifndef NC
  #define NC 3
#endif

// default number of spins
#ifndef NS
  #define NS 4
#endif

// number of dimensions
#ifndef ND
  #define ND 4
#endif

#define NSNS ( NS*NS )

#define NCNC ( NC*NC )

#define PLAQ_AND_TRACE_TOL 1E-6

#define PREC_TOL 1E-14

#endif
