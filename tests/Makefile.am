## simple makefile for the unit tests
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = MESONS
MESONS_SOURCES = mesonfile.c

## compile against the lib installed somewhere
MESONS_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/HEADERS/
MESONS_LDADD = ${TOPDIR}/src/libCORR.a ${LDFLAGS}
