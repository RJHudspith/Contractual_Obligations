/**
   @mainpage The CORR library

   @author Anthony Francis
   @author Renwick Hudspith

   @section Introduction

   Fermion propagator contraction code. We required a code that would take output propagator solutions from different codes (sometimes chiral and sometimes non-relativistic) and contracted these to form the Mesons and Baryons that we all love.

   As the code has progressed we have added some features, here is a (non-extensive) list of features,

   > Mesons with non zero momentum projection

   > Baryons ( flavour diagonal and non ) with non-zero momentum projection

   > Conserved-Local and Local-Local Wilson currents

   > Weak Matrix Elements

   This code has been written to work for Wall and point sources. Some more technical details are:

   > SSE versions of lowest level linear algebra

   > FFTW incorporation for all Fourier Transforms

   > Openmp threading

   > Full c99 compliance and clean compilation under -wpedantic with clang and gcc
*/
