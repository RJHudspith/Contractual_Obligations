## makes the documentation
doxygen:
	doxygen Doxyfile

.PHONY: 
	doxygen

documentation:
	pdflatex intro.tex && mv -v *.pdf ${DOCPREFIX}/docs/ && rm *.out *.log *.aux *~

.PHONY:
	documentation