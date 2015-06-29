## simple makefile for the unit tests
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = BARYONS MESONS MOMAVG UNIT VPFREAD

## meson reader
BARYONS_SOURCES = baryonfile.c reader.c
BARYONS_CFLAGS = -I${TOPDIR}/src/HEADERS/ 
BARYONS_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## meson reader
MESONS_SOURCES = mesonfile.c reader.c
MESONS_CFLAGS = -I${TOPDIR}/src/HEADERS/
MESONS_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## momentum average reader
MOMAVG_SOURCES = momavg.c
MOMAVG_CFLAGS = -I${TOPDIR}/src/HEADERS/
MOMAVG_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## (non-exhaustive) unit test compilation
UNIT_SOURCES = unit.c contract_tests.c matops_tests.c \
	spinmatrix_tests.c spinor_tests.c
UNIT_CFLAGS = -I${TOPDIR}/src/HEADERS/
UNIT_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS} 

## vpf reader
VPFREAD_SOURCES = vpffile.c
VPFREAD_CFLAGS = -I${TOPDIR}/src/HEADERS/
VPFREAD_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS} 