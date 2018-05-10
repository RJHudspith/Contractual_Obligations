CO
==

Thread-parallel Fermion propagator contraction code written 
by A. Francis and J. Hudspith.

INSTALLATION
============

You will need to run these on the command line to get the configure script

aclocal && autoreconf -ivf

Then it is the usual automake stuff

./configure --prefix={prefdir} --with-fftw={path to fftw}

The code automatically compiles in ND = NS = 4, NC=3. To use something
different you would need to call

--with=NC= --with-ND= --with-NS=

Some codes (such as the Baryons/Pentas are not yet implemented to be NC-generic)

USAGE
=====

./CORR -i infile # -g {gauge file}

(the stuff after the # is optional, used for the VPF)

The infile is then used to choose the things you wish to contract

An auxillary binary ./MESONS can be used to inspect various meson files

./MESON {file name} GSRC1,GSNK1 GSRC2,GSNK2 ....

INTRINSICS
==========

By default, if configure can find SSE headers (immintrin.h) it will
switch on the code that uses the vector intrinsics. We align to 16 bytes, and have access to FMA and some small AVX2 usage.

FUNCTIONALITY
=============

As of writing (25/06/15) the code supports

Meson contractions at zero and non-zero spatial momentum
Conserved-local and local-local Wilson currents (as long as the
following gauge field configuration format is used: NERSC, HiREP or LIME)
Baryon contractions for SU(2) and SU(3)
Limited static quark propagators
!Untested! Weak matrix element contractions

We now have Tetraquarks and Pentaquarks

TESTS
=====

We have unit test coverage on our linear algebra, testing the most basic
matrix operations for correctness. So all spinor operations and all
color matrix operations and all simple contractions are
covered.

DOCUMENTATION
=============

Our code is doxygen documented. If available, the command

make doxygen

will make it all in the library directory (in ${Installation dir})
and will create the doxygen webpage in

${Installation dir}/docs/html/index.html

and can be opened with your internet browser.

the command

make documentation

will attempt to compile (using pdflatex) our documentation too, this
will create the file 

${Installation dir}/docs/Corr.pdf