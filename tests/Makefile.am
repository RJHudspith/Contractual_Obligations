## simple makefile for the unit tests
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

## binaries generated
bin_PROGRAMS = UNIT

## (non-exhaustive) unit test compilation
UNIT_SOURCES = unit.c contract_tests.c matops_tests.c \
	spinmatrix_tests.c spinor_tests.c \
	bar_projections_tests.c bar_ops_tests.c \
	halfspinor_tests.c \
	tetra_contractions_tests.c \
	gamma_tests.c utils_tests.c \
	SSE_tests.c
UNIT_CFLAGS = -I${TOPDIR}/src/HEADERS/
UNIT_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}
