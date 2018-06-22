TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libCORR.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

## c files in ./BARYONS/
BARYONFILES=./BARYONS/bar_contractions.c ./BARYONS/bar_projections.c \
	./BARYONS/bar_ops.c ./BARYONS/bar_ops_SSE.c \
	./BARYONS/baryons_uuu.c ./BARYONS/baryons_uud.c \
	./BARYONS/baryons_uds.c ./BARYONS/wrap_baryons.c

## c files in ./DIQUARKS/
DIQUARKFILES=./DIQUARKS/diquark.c ./DIQUARKS/diquark_degen.c \
	./DIQUARKS/wrap_diquarks.c ./DIQUARKS/diquark_contraction.c

## c files in ./GEOM/
GEOMFILES=./GEOM/gammas.c ./GEOM/geometry.c ./GEOM/plan_ffts.c

## c files in ./FIELD/
FIELDFILES=./FIELD/plaqs_links.c

## c files in ./IO/
IOFILES=./IO/CERN.c ./IO/chklat_stuff.c ./IO/correlators.c ./IO/cut_output.c \
	./IO/HIREP.c ./IO/io.c ./IO/input_reader.c \
	./IO/input_baryons.c ./IO/input_general.c ./IO/input_mesons.c \
	./IO/input_pentas.c ./IO/input_tetras.c ./IO/input_VPF.c \
	./IO/input_WME.c ./IO/read_config.c ./IO/read_headers.c \
	./IO/read_propheader.c \
	./IO/readers.c ./IO/Scidac.c ./IO/XML_info.c

## c files in ./LINALG/
LINALGFILES= \
	./LINALG/contractions.c ./LINALG/contractions_SSE.c \
	./LINALG/halfspinor_ops.c \
	./LINALG/matrix_ops.c ./LINALG/matrix_ops_SSE.c \
	./LINALG/mmul.c ./LINALG/mmul_SSE.c \
	./LINALG/Ospinor.c \
	./LINALG/spinor_ops.c ./LINALG/spinor_ops_SSE.c \
	./LINALG/spinmatrix_ops.c ./LINALG/spinmatrix_ops_SSE.c

## c files in ./MEAS/
MEASFILES=./MEAS/mesons.c ./MEAS/mesons_offdiag.c ./MEAS/wrap_mesons.c

## c files in ./NRQCD
NRQCDFILES=./NRQCD/nrqcd.c ./NRQCD/sources.c ./NRQCD/spin_independent.c \
	./NRQCD/clover.c ./NRQCD/grad.c ./NRQCD/grad_2.c ./NRQCD/grad_4.c \
	./NRQCD/spin_dependent.c ./NRQCD/evolve.c

## c files in ./PENTA/
PENTAFILES=./PENTA/contract_O1O1.c ./PENTA/contract_O1O2.c \
	./PENTA/contract_O1O3.c ./PENTA/contract_O2O1.c \
	./PENTA/contract_O2O2.c ./PENTA/contract_O2O3.c \
	./PENTA/contract_O3O1.c ./PENTA/contract_O3O2.c \
	./PENTA/contract_O3O3.c \
	./PENTA/penta_contractions.c \
	./PENTA/penta_udusb.c ./PENTA/wrap_pentas.c

## c files in ./TETRA/
TETRAFILES=./TETRA/tetra_usbb.c ./TETRA/wrap_tetras.c \
	./TETRA/tetra_contractions.c ./TETRA/tetra_udbb.c \
	./TETRA/tetra_udcb.c ./TETRA/tetra_uscb.c

## c files in ./UTILS/
UTILSFILES=./UTILS/basis_conversions.c ./UTILS/corr_malloc.c \
	./UTILS/corr_sort.c ./UTILS/crc32.c ./UTILS/crc32c.c \
	./UTILS/cut_routines.c ./UTILS/GLU_bswap.c \
	./UTILS/GLU_timer.c ./UTILS/gramschmidt.c \
	./UTILS/progress_bar.c ./UTILS/setup.c \
	./UTILS/par_MWC_4096.c ./UTILS/par_rng.c

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
	${DIQUARKFILES} \
	${GEOMFILES} \
	${FIELDFILES} \
	${IOFILES} \
	${LINALGFILES} \
	${MEASFILES} \
	${NRQCDFILES} \
	${PENTAFILES} \
	${TETRAFILES} \
	${UTILSFILES} \
	${VPFFILES} \
	${WMEFILES}

## pattern substitute in the headers
include_HEADERS = \
	$(BARYONFILES:./BARYONS/%.c=./HEADERS/%.h) \
	$(DIQUARKFILES:./DIQUARKS/%.c=./HEADERS/%.h) \
	$(GEOMFILES:./GEOM/%.c=./HEADERS/%.h) \
	$(FIELDFILES:./FIELD/%.c=./HEADERS/%.h) \
	$(IOFILES:./IO/%.c=./HEADERS/%.h) \
	$(LINALGFILES:./LINALG/%.c=./HEADERS/%.h) \
	$(MEASFILES:./MEAS/%.c=./HEADERS/%.h) \
	$(NRQCDFILES:./NRQCD/%.c=./HEADERS/%.h) \
	$(PENTAFILES:./PENTA/%.c=./HEADERS/%.h) \
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
CORR_LDADD = libCORR.a ${LDFLAGS} -lm

endif
