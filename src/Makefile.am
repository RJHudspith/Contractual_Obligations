TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libCORR.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

## c files in ./BARYONS/
BARYONFILES= ./BARYONS/baryons.c ./BARYONS/baryons2.c ./BARYONS/baryons3.c \
	./BARYONS/wrap_baryons.c

## c files in ./GEOM/
GEOMFILES=./GEOM/gammas.c ./GEOM/geometry.c ./GEOM/plan_ffts.c

## c files in ./IO/
IOFILES=./IO/chklat_stuff.c ./IO/cut_output.c ./IO/io.c ./IO/input_reader.c \
	./IO/input_baryons.c ./IO/input_general.c ./IO/input_mesons.c \
	./IO/input_tetras.c ./IO/input_VPF.c ./IO/input_WME.c \
	./IO/read_config.c ./IO/read_headers.c ./IO/read_propheader.c \
	./IO/readers.c ./IO/Scidac.c ./IO/XML_info.c

## c files in ./MEAS/
MEASFILES=./MEAS/correlators.c ./MEAS/mesons.c ./MEAS/mesons_offdiag.c \
	./MEAS/plaqs_links.c ./MEAS/wrap_mesons.c

## c files in ./TETRA/
TETRAFILES=./TETRA/tetraquark.c ./TETRA/wrap_tetras.c

## c files in ./UTILS/
UTILSFILES=./UTILS/bar_contractions.c ./UTILS/bar_ops.c ./UTILS/bar_ops_SSE.c \
	./UTILS/basis_conversions.c ./UTILS/contractions.c \
	./UTILS/contractions_SSE.c ./UTILS/corr_malloc.c ./UTILS/corr_sort.c \
	./UTILS/crc32.c ./UTILS/cut_routines.c \
	./UTILS/GLU_bswap.c ./UTILS/GLU_timer.c ./UTILS/gramschmidt.c \
	./UTILS/matrix_ops.c ./UTILS/matrix_ops_SSE.c \
	./UTILS/spinor_ops.c ./UTILS/spinor_ops_SSE.c

## c files in ./VPF
VPFFILES=./VPF/cl_diagonal.c ./VPF/cl_offdiagonal.c ./VPF/currents.c \
	./VPF/ll_diagonal.c ./VPF/ll_offdiagonal.c \
	./VPF/momspace_PImunu.c ./VPF/PImunu_projections.c \
	./VPF/tmoments_PImunu.c ./VPF/WardIdentity.c ./VPF/wrap_VPF.c

## c files in ./WME/
WMEFILES=./WME/WME.c ./WME/wrap_WME.c

## c file sources for the library 
libCORR_a_SOURCES = \
	${BARYONFILES} \
	${GEOMFILES} \
	${IOFILES} \
	${MEASFILES} \
	${TETRAFILES} \
	${UTILSFILES} \
	${VPFFILES} \
	${WMEFILES}

## pattern substitute in the headers
include_HEADERS = \
	$(BARYONFILES:./BARYONS/%.c=./HEADERS/%.h) \
	$(GEOMFILES:./GEOMETRY/%.c=./HEADERS/%.h) \
	$(IOFILES:./IO/%.c=./HEADERS/%.h) \
	$(MEASFILES:./MEAS/%.c=./HEADERS/%.h) \
	$(TETRAFILES:./TETRA/%.c=./HEADERS/%.h) \
	$(UTILSSFILES:./UTILS/%.c=./HEADERS/%.h) \
	$(VPFFILES:./VPF/%.c=./HEADERS/%.h) \
	$(WMEFILES:./WME/%.c=./HEADERS/%.h) \
	./HEADERS/common.h \
	../config.h 

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = CORR

CORR_SOURCES = corr.c
CORR_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
CORR_LDADD = libCORR.a ${LDFLAGS}

endif