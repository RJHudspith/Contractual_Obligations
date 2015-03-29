/**
   @file common.h
   @brief some common headers and some defines, this will be changed
 */
#ifndef COMMON_H
#define COMMON_H

#include "../config.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// SSE2 header
#if (defined HAVE_EMMINTRIN_H)
  #include <emmintrin.h>
#endif

// SSE3 header
#if (defined HAVE_PMMINTRIN_H)
  #include <pmmintrin.h>
#endif

// fftw instructions
#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#endif

// 
#ifndef WORDS_BIGENDIAN
  #define WORDS_BIGENDIAN 0
#endif

#include "definitions.h"

#include "enum.h"

#include "structs.h"

// important info gets stored here
extern struct latt_info Latt ;

#endif
