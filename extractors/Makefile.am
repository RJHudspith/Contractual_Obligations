## simple makefile for extractors
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS=BARYONS CAVG MESONS MOMAVG VPFREAD TETAVG PROPOPS PENTAVG

## baryon reader
BARYONS_SOURCES=baryonfile.c reader.c
BARYONS_CFLAGS=-I${TOPDIR}/src/HEADERS/ 
BARYONS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## baryon reader
CAVG_SOURCES=cavg.c reader.c
CAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/ 
CAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## meson reader
MESONS_SOURCES=mesonfile.c reader.c
MESONS_CFLAGS=-I${TOPDIR}/src/HEADERS/
MESONS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## momentum average reader
MOMAVG_SOURCES=momavg.c reader.c
MOMAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/
MOMAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## vpf reader
VPFREAD_SOURCES=vpffile.c
VPFREAD_CFLAGS=-I${TOPDIR}/src/HEADERS/
VPFREAD_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## tet averager
TETAVG_SOURCES=tet_avg.c reader.c
TETAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/
TETAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## tet averager
PENTAVG_SOURCES=pent_avg.c reader.c
PENTAVG_CFLAGS=-I${TOPDIR}/src/HEADERS/
PENTAVG_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm

## prop operators
PROPOPS_SOURCES=PropOps.c
PROPOPS_CFLAGS=-I${TOPDIR}/src/HEADERS/
PROPOPS_LDADD=${TOPDIR}/src/libCORR.a ${LDFLAGS} -lm
