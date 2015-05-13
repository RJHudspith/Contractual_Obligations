## simple makefile for the unit tests
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = MESONS UNIT VPFREAD

## meson reader
MESONS_SOURCES = mesonfile.c
MESONS_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
MESONS_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}

## (non-exhaustive) unit test compilation
UNIT_SOURCES = unit.c contract_tests.c matops_tests.c spinor_tests.c
UNIT_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
UNIT_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}

## vpf reader
VPFREAD_SOURCES = vpffile.c
VPFREAD_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
VPFREAD_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}