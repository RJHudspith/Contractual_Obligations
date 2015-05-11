## makes the documentation
doxygen:
	doxygen Doxyfile

.PHONY: 
	doxygen

documentation:
	pdflatex Corr.tex && bibtex Corr.aux && pdflatex Corr.tex && pdflatex Corr.tex && mv -v *.pdf ${DOCPREFIX}/docs/ && rm *.out *.log *.aux *~ *.bbl *.blg

.PHONY:
	documentation