TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libCORR.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

## c files in ./BARYONS/
BARYONFILES=./BARYONS/wrap_baryons.c

## c files in ./GEOM/
GEOMFILES=./GEOM/gammas.c ./GEOM/geometry.c ./GEOM/plan_ffts.c

## c files in ./IO/
IOFILES=./IO/chklat_stuff.c ./IO/io.c ./IO/input_reader.c ./IO/read_config.c \
	./IO/readers.c ./IO/read_headers.c ./IO/Scidac.c ./IO/XML_info.c \
	./IO/read_propheader.c ./IO/cut_output.c

## c files in ./MEAS/
MEASFILES=./MEAS/correlators.c ./MEAS/plaqs_links.c ./MEAS/brutal_mesons.c \
	./MEAS/wrap_mesons.c ./MEAS/mesons.c

## c files in ./VPF
VPFFILES=./VPF/wrap_VPF.c ./VPF/PImunu_projections.c ./VPF/currents.c \
	./VPF/WardIdentity.c ./VPF/local_local.c ./VPF/conserved_local.c \
	./VPF/momspace_PImunu.c ./VPF/tmoments_PImunu.c

## c files in ./UTILS/
UTILSFILES=./UTILS/GLU_timer.c ./UTILS/gramschmidt.c ./UTILS/matrix_ops.c \
	./UTILS/spinor_ops.c ./UTILS/basis_conversions.c ./UTILS/contractions.c \
	./UTILS/cut_routines.c ./UTILS/crc32.c ./UTILS/GLU_bswap.c

## c files in ./WME/
WMEFILES=./WME/WME.c ./WME/wrap_WME.c

## now I think I can probably pattern substitute all the c files into headers TODO

include_HEADERS = \
	../config.h ./HEADERS/chklat_stuff.h ./HEADERS/common.h ./HEADERS/correlators.h \
	./HEADERS/conserved_local.h ./HEADERS/crc32.h ./HEADERS/definitions.h \
	./HEADERS/enum.h ./HEADERS/gammas.h ./HEADERS/mesons.h \
	./HEADERS/geometry.h ./HEADERS/GLU_bswap.h ./HEADERS/GLU_timer.h \
	./HEADERS/gramschmidt.h ./HEADERS/input_reader.h ./HEADERS/io.h \
	./HEADERS/matrix_ops.h	./HEADERS/momspace_PImunu.h \
	./HEADERS/plaqs_links.h ./HEADERS/read_config.h ./HEADERS/readers.h \
	./HEADERS/read_headers.h ./HEADERS/Scidac.h ./HEADERS/structs.h ./HEADERS/wall_mesons.h \
	./HEADERS/XML_info.h ./HEADERS/read_propheader.h ./HEADERS/brutal_mesons.h \
	./HEADERS/currents.h ./HEADERS/spinor_ops.h ./HEADERS/wrap_mesons.h \
	./HEADERS/basis_conversions.h ./HEADERS/contractions.h ./HEADERS/WardIdentity.h \
	./HEADERS/tmoments_PImunu.h ./HEADERS/cut_routines.h ./HEADERS/cut_output.h \
	./HEADERS/wrap_VPF.h ./HEADERS/local_local.h ./HEADERS/PImunu_projections.h \
	./HEADERS/WME.h ./HEADERS/wrap_WME.h ./HEADERS/wrap_baryons.h 

libCORR_a_SOURCES = \
	${BARYONFILES} \
	${GEOMFILES} \
	${IOFILES} \
	${MEASFILES} \
	${UTILSFILES} \
	${VPFFILES} \
	${WMEFILES}

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = CORR

CORR_SOURCES = corr.c
CORR_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
CORR_LDADD = libCORR.a ${LDFLAGS}

endif