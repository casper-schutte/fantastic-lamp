all: main.pdf


main.pdf: main.tex document.bib *.tex
	pdflatex main </dev/null
	bibtex main </dev/null
	pdflatex main.tex </dev/null
	pdflatex main.tex </dev/null

clean:
	rm -f main.aux main.toc main.pdf main.log main.bbl main.blg
