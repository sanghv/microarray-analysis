#
# ACM requires ptmb8t fonts. Under aptitude package manager in debian systems:
#	sudo apt-get install texlive-fonts-recommended 
#
#	Created: cw
#
##################################################################################

filename=sigproc-sp

# Generate only pdf
pdf: ${filename}.tex 
	pdflatex ${filename}.tex
	bibtex ${filename}.aux || true
	pdflatex ${filename}.tex
	pdflatex ${filename}.tex

# Read the pdf file.
aread:
	evince ${filename}.pdf &

# get rid of all generated files
clobber:
	rm -f ${filename}.ps ${filename}.pdf ${filename}.log ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg ${filename}.dvi.pdf

# keep processed documents. Get rid of everything else.
clean:
	${filename}.log ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg ${filename}.dvi.pdf


