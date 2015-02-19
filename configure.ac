## GLU Automake instructions outline ##

## --prefix=  :: tell it where to install the binary
## --with-NC= :: number of colors, if not specified defaults to 3
## --with-ND= :: number of dimensions, if not specified defaults to 4

AC_INIT([CORR], [0.1], [renwick.james.hudspith@gmail.com,afrancis.heplat@googlemail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
##AM_PROG_AR
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AC_C_BIGENDIAN 

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
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
