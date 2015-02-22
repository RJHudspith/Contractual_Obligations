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

#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#endif

#define WORDS_BIGENDIAN 0

#include "definitions.h"

#include "enum.h"

#include "structs.h"

// important info gets stored here
extern struct latt_info Latt ;

#endif
