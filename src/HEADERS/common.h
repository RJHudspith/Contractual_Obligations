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

// I need to think about the logic here
#if (defined HAVE_IMMINTRIN_H)
  #include <immintrin.h>
  #if (defined __SSE2__ ) || (defined __SSE3__ )
    #define HAVE_EMMINTRIN_H
    #include "SSE2_OPS.h"
  #endif
#endif

// needed a wrapper for undefined posix_memalign
#include "corr_malloc.h"

// fftw instructions
#ifdef HAVE_FFTW3_H
  #include <fftw3.h>
#endif

// wrap these openmp functions
#if (defined _OPENMP ) && (defined HAVE_OMP_H )
  #include <omp.h>
  #define get_CORR_thread() omp_get_thread_num()
#else
  #define get_CORR_thread() (0)
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

// gauge field is global !!
extern struct site *lat ;

#endif
