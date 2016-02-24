## simple makefile for extractors
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS=BARYONS MESONS MOMAVG VPFREAD

## meson reader
BARYONS_SOURCES=baryonfile.c reader.c
BARYONS_CFLAGS=-I${TOPDIR}/src/HEADERS/ 
BARYONS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## meson reader
MESONS_SOURCES=mesonfile.c reader.c
MESONS_CFLAGS=-I${TOPDIR}/src/HEADERS/
MESONS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## momentum average reader
MOMAVG_SOURCES=momavg.c reader.c
MOMAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/
MOMAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## vpf reader
VPFREAD_SOURCES=vpffile.c
VPFREAD_CFLAGS=-I${TOPDIR}/src/HEADERS/
VPFREAD_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} 