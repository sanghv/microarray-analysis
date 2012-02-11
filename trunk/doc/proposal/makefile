#
#	Basic make file:
#		make  		 --- builds
#		make clean   --- cleans stuff up
#		make aread	 --- Read the pdf 
#
################################################3


filename=proposal

pdf: ${filename}.tex *.bib
	pdflatex ${filename}.tex
	bibtex ${filename}.aux
	pdflatex ${filename}.tex
	pdflatex ${filename}.tex

aread:
	evince ${filename}.pdf &

clean:
	rm -f ${filename}.ps ${filename}.pdf ${filename}.log ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg ${filename}.dvi.pdf
