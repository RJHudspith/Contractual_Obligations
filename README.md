CO
==

Fermion propagator contraction code

INSTALLATION
============

You will need to run these on the command line to get the configure script

aclocal && autoreconf -ivf

Then it is the usual automake stuff

./configure --prefix={prefdir} --with-fftw={path to fftw}

USAGE
=====

./CORR -i infile # -g {gauge file}

(the stuff after the # is optional, used for the VPF)

The infile is then used to choose the things you wish to contract

An auxillary binary ./MESONS can be used to inspect various meson files

./MESON {file name} GSRC1,GSNK1 GSRC2,GSNK2 ....