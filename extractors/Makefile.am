## simple makefile for extractors
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS=BARYONS MESONS MOMAVG VPFREAD TETAVG PROPOPS

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

## tet averager
TETAVG_SOURCES=tet_avg.c reader.c
TETAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/
TETAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS}

## prop operators
PROPOPS_SOURCES=PropOps.c
PROPOPS_CFLAGS=-I${TOPDIR}/src/HEADERS/
PROPOPS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS}
