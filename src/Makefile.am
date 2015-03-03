TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libCORR.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

include_HEADERS = \
	../config.h ./HEADERS/chklat_stuff.h ./HEADERS/common.h ./HEADERS/correlators.h \
	./HEADERS/conserved_local.h \
	./HEADERS/crc32.h ./HEADERS/definitions.h ./HEADERS/enum.h ./HEADERS/gammas.h \
	./HEADERS/geometry.h ./HEADERS/GLU_bswap.h ./HEADERS/GLU_timer.h \
	./HEADERS/gramschmidt.h ./HEADERS/input_reader.h ./HEADERS/io.h \
	./HEADERS/matrix_ops.h	./HEADERS/momspace_PImunu.h ./HEADERS/mesons.h \
	./HEADERS/plaqs_links.h ./HEADERS/read_config.h ./HEADERS/readers.h \
	./HEADERS/read_headers.h ./HEADERS/Scidac.h ./HEADERS/structs.h ./HEADERS/wall_mesons.h \
	./HEADERS/XML_info.h ./HEADERS/read_propheader.h ./HEADERS/brutal_mesons.h \
	./HEADERS/currents.h ./HEADERS/spinor_ops.h

libCORR_a_SOURCES = \
	./GEOM/gammas.c ./GEOM/geometry.c ./GEOM/plan_ffts.c \
	./IO/chklat_stuff.c ./IO/io.c ./IO/input_reader.c ./IO/read_config.c \
	./IO/readers.c ./IO/read_headers.c ./IO/Scidac.c ./IO/XML_info.c \
	./IO/read_propheader.c \
	./MEAS/correlators.c ./MEAS/mesons.c ./MEAS/wall_mesons.c \
	./MEAS/conserved_local.c ./MEAS/momspace_PImunu.c \
	./MEAS/plaqs_links.c ./UTILS/crc32.c ./UTILS/GLU_bswap.c \
	./MEAS/brutal_mesons.c ./MEAS/currents.c \
	./UTILS/GLU_timer.c ./UTILS/gramschmidt.c ./UTILS/matrix_ops.c \
	./UTILS/spinor_ops.c

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = CORR

CORR_SOURCES = corr.c
CORR_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
CORR_LDADD = libCORR.a ${LDFLAGS}

endif