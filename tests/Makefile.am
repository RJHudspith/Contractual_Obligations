## simple makefile for the unit tests
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = MESONS UNIT
MESONS_SOURCES = mesonfile.c
UNIT_SOURCES = unit.c contract_tests.c matops_tests.c spinor_tests.c

## compile against the lib installed somewhere
MESONS_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
MESONS_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}

## (non-exhaustive) unit test compilation
UNIT_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
UNIT_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}