## GLU Automake instructions outline ##

## --prefix=  :: tell it where to install the binary
## --with-NC= :: number of colors, if not specified defaults to 3
## --with-ND= :: number of dimensions, if not specified defaults to 4

AC_INIT([CORR], [0.2], [renwick.james.hudspith@gmail.com,afrancis.heplat@googlemail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
m4_pattern_allow([AM_PROG_AR])
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AC_C_BIGENDIAN
AM_PROG_AR

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

## prefix!!
AC_DEFINE_UNQUOTED([HAVE_PREFIX],["$prefix"],[Look for where we prefix to])
	  
## Tells us where we are prefixed to
AC_MSG_NOTICE([Configuring for location "$prefix"])

## If we prefix away from the default /usr/ then create binary GLU
## will still try and create shared object libraries there
AM_CONDITIONAL([PREF],[ test "$prefix" == "NONE" ])
if test "$prefix" == "NONE" ; 
then 
    AC_MSG_NOTICE([Not creating binary "CORR"])
    AM_CONDITIONAL([DOX],[false])
    AM_CONDITIONAL([DOT],[false])
    AM_CONDITIONAL([PDFLTX],[false])
else
    AC_MSG_NOTICE([Creating binary "CORR" in $prefix/bin])
    ## check for doxygen
    AC_CHECK_PROG([DOXYGEN],[doxygen],["yes"],["NONE"])
    AM_CONDITIONAL([DOX],[ test "$DOXYGEN" == "NONE" ])
    ## check for dot
    AC_CHECK_PROG([GRAPHVIZ],[dot],["yes"],["NONE"])
    AM_CONDITIONAL([DOT],[ test "$GRAPHVIZ" == "NONE" ])
    ## check for pdflatex
    AC_CHECK_PROG([PDFLATEX],[pdflatex],["yes"],["NONE"])
    AM_CONDITIONAL([PDFLTX],[ test "$PDFLATEX" == "NONE" ])
    mkdir $prefix  
    mkdir $prefix/bin
    mkdir $prefix/docs/ ## new documentation folder ...
    ## copy the input file and a couple of small tests to bin	
    cp -v ./CONF/* $prefix/bin/
fi
AC_SUBST(PREF)

## Compile for a set NC, default is NC=3
AC_ARG_WITH([NC],
	AC_HELP_STRING([--with-NC=],
		       [Compile in the gauge group number NC]),
            [
	    AC_MSG_NOTICE([User specified NC, compiling for SU(${with_NC})])
	    AC_DEFINE_UNQUOTED([NC], [${with_NC}] , [Compiled for SU(NC)] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified NC, default to SU(3)])	
	    ])

## Compile for a set ND, default is ND=4
AC_ARG_WITH([ND],
	AC_HELP_STRING([--with-ND=],
		       [Compile in the number of dimensions ND]),
            [
	    AC_MSG_NOTICE([User specified ND, compiling for ND=${with_ND}])
	    AC_DEFINE_UNQUOTED([ND], [${with_ND}] , [Compiled for ND] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified ND, default to 4])	
	    ])

## Compile for a set NS, default is ND=4
AC_ARG_WITH([NS],
	AC_HELP_STRING([--with-NS=],
		       [Compile in the number of spins NS]),
            [
	    AC_MSG_NOTICE([User specified NS, compiling for NS=${with_NS}])
	    AC_DEFINE_UNQUOTED([NS], [${with_NS}] , [Compiled for NS] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified NS, default to 4])	
	    ])

## Finally we have routines for this node or many machines
## historically targeted at "CONDOR", default is CONDOR_MODE
AC_ARG_ENABLE([notcondor],
	[ --enable-notcondor Allows for the saving of FFT plans and some simple hashing and whatever],
	[case "${enableval}" in
	yes) ## Check the condor enable creates the dirs $prefix/Local/(Wisdom && Moments)
	     ## used for some light caching of previous results
             if [ test "$prefix" == "NONE" ] ;then
             	AC_MSG_NOTICE([Cannot save wisdom as we are not prefixed, defining CONDOR_MODE by default])
             else
		AC_MSG_NOTICE([FFTW's WISDOM routines being used, targeted specifically for this node])
		AC_DEFINE([NOT_CONDOR_MODE],,[targeted, has been seen to be quicker but do not distribute the binaries!])
		mkdir $prefix/Local/	
		mkdir $prefix/Local/Wisdom
		mkdir $prefix/Local/Moments 
             fi
	     AC_MSG_NOTICE([FFTW plans created ONCE. Subsequent transforms are read from Local/Wisdom/])
	     ;;
	no) AC_MSG_NOTICE([default condor routines called, i.e. FFTW plans created here and now]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-notcondor]) ;;
	esac],[])

## printfs debug information
AC_ARG_ENABLE([debug],
	[ --enable-debug prints to stdout some general debugging information],
	[case "${enableval}" in
	yes) AC_DEFINE([DEBUG],,[general debug information])
	     AC_MSG_NOTICE([Printing debug information to stdout]) ;;	
	no) AC_MSG_NOTICE([Not printing debug information to stdout]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
	esac],[])

## TAKEN LIKE MANY OTHER THINGS IN THIS CODE FROM MY LIBRARY GLU - J
######################################### FFTW3 ###############################################
##                                                                                           ##
##    This one is tricky, looks first for the static link to FFTW library, if not found      ##
##    then the shared. If --with-fftw= then it looks in /usr/lib for the shared              ##
##    otherwise we miss out on all of the Fourier transform stuff like Fourier acceleration, ##
##    Momentum space gluons and quenched U(1). Not a necessity.                              ##
##                                                                                           ##
###############################################################################################

## I have now written the code to utilise OpenMP parallel FFT routines
AC_ARG_ENABLE([OMP_FFTW],
	[ --enable-OMP_FFTW Utilise FFTW's OpenMP parallel FFTs],
	[case "${enableval}" in
	yes) OMP_FFTW=true ;;
	no) OMP_FFTW=false ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-OMP_FFTW]) ;;
	esac],[OMP_FFTW=false])

## Look for openmp, cannot turn on the omp fftw routines without it!
AC_CHECK_HEADERS([omp.h],[],[OMP_FFTW=false])

## Look for fftw in the useer specified place
fftw=false 
AC_ARG_WITH([fftw],
	AC_HELP_STRING([--with-fftw=<name>],
		       [Specify FFTW location for the GF]),
            [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])
	    
	    FFTW=${with_fftw}

	    ## are we compiling for single precision or not?
	    if test x$single = xtrue ; then
	       float_version=f
	    else	
	       float_version=""
	    fi

	    AC_CHECK_HEADERS([${FFTW}include/fftw3.h],[
	    AC_MSG_NOTICE([User defined fftw library @ $FFTW])
	    AC_MSG_NOTICE([Attempting linking to the user's FFTW lib ${FFTW}lib/libfftw3${float_version}])
	
	    ## if we use the user one, only link to the static lib
	    if [ test -a ${FFTW}lib/libfftw3${float_version}.a ]; then
		AC_MSG_NOTICE([Statically linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include"
	    	LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}.a"
		if test x$OMP_FFTW = xtrue ; then
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.a ; then
		      LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}_omp.a"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    ## Have a look for the shared object library ...
            elif [ test -a ${FFTW}lib/libfftw3${float_version}.so ]; then
		AC_MSG_NOTICE([Static link not found !! Shared object (.so) linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include -L${FFTW}lib/"
	    	LDFLAGS="$LDFLAGS -lfftw3${float_version}"
		if test x$OMP_FFTW = xtrue ; then
		   ## check to see if it is there otherwise don't include it
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.so ; then
		      LDFLAGS="$LDFLAGS -lfftw3${float_version}_omp"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
            else
	    	AC_MSG_ERROR([fftw3.h found, but cannot find the shared or static libraries!])
	    fi
	    ]	    
	    , [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])

	    ## get desperate here and check for the header in the usual places
	    ## as omp routines are special I doubt they will be installed, maybe should do a special check
	    if test x$local_fftw = "xtrue" ; then
	       LDFLAGS="$LDFLAGS -lfftw3${float_version}"
	       AC_MSG_NOTICE([Using the system's fftw])	       
	       AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    else
	       AC_MSG_NOTICE([No local fftw3.h found either!])
	    fi
     ]) ## closes the check headers ...
]
## closes with-FFTW
)
AM_CONDITIONAL([DEFFFTW],[$fftw])
AC_SUBST(FFTW)

## Compile a tetra block size
AC_ARG_WITH([TETRA_NBLOCK],
	AC_HELP_STRING([--with-TETRA_NBLOCK=],
		       [Block matrix tetra contraction size]),
            [
	    AC_MSG_NOTICE([User specified TETRA_NBLOCK, compiling for (${with_TETRA_NBLOCK})])
	    AC_DEFINE_UNQUOTED([TETRA_NBLOCK], [${with_TETRA_NBLOCK}] , [Compiled TETRA_NBLOCK] )	
	    ],[
	    AC_MSG_NOTICE([User unspecified TETRA_NBLOCK, code defaults to 1])	
])

## C standard libraries
AC_CHECK_HEADERS([math.h],[LDFLAGS="$LDFLAGS -lm" ])

## C standard library 
AC_CHECK_HEADERS([stdint.h],[])

## check for the unix memory stuff,
AC_CHECK_HEADERS([unistd.h], 
	[AC_MSG_NOTICE([We are attempting to be memory-safe in our operation.])],[])

## Check for sys/time.h
AC_CHECK_HEADERS([sys/time.h],
	[AC_MSG_NOTICE([We are intrinsically timing our operations.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([time.h],
	[AC_MSG_NOTICE([We are providing an accurate date.])],[])

## SSE headers
AC_CHECK_HEADERS([immintrin.h],[],[])

## Some compiler checks
## My code is littered with consts
AC_C_CONST([])

## Also use the restrict flag a lot
AC_C_RESTRICT([])

## Many of my functions are "inlined" {compiler's choice so maybe not} ,
## hence the large binary
AC_C_INLINE([])

## Set up config.h, one of the most important headers of them all
AC_CONFIG_HEADERS([config.h])

## where are my files?
AC_CONFIG_FILES([
Makefile
src/Makefile
extractors/Makefile
tests/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
