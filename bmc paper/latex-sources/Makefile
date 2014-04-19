.SILENT:

all: bmc_article.pdf

bmc_article.bbl: bmc_article.bib
	latex bmc_article.tex
	bibtex bmc_article

bmc_article.pdf: $(wildcard *.tex) $(wildcard img/*.eps) $(wildcard img/*.pdf) bmc_article.bbl
	latex bmc_article.tex
	latex bmc_article.tex
	dvips bmc_article.dvi
	ps2pdf bmc_article.ps

clean:
	rm -f *.ps *.pdf *.dvi *.bbl *.blg *.aux *.log *.toc *.backup *.fdb_pdflatexmk \#*  *~
