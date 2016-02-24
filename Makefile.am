## everything is in src and tests
SUBDIRS = src tests extractors

ACLOCAL_AMFLAGS = -I m4

####################################### AUTOMATIC DOCUMENTATION GENERATION #####################################

## Generates the doxygen for the code
if !DOX
## if we have "dot" we generate call graphs
if !DOT
  have_dot=YES
else
  havedot=NO
endif
doxygen:
	sed s~MYOUTDIR~"${prefix}/docs"~g ./docs/Doxyfile.prov | sed s~MYDOT~${have_dot}~g > ./docs/Doxyfile 
	(cd docs && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1

.PHONY: 
	doxygen
endif

## if we have pdflatex, I allow for the rule to autogenerate the documentation
## from the latex source
if !PDFLTX
documentation:
	(cd docs && $(MAKE) DOCPREFIX="${prefix}" $(AM_MAKEFLAGS) $@) || exit 1
.PHONY: 
	documentation
endif