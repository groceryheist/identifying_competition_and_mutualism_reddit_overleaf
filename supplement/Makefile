#!/usr/bin/make
all: $(patsubst %.Rtex,%.pdf,$(wildcard *.Rtex)) 

ecological_models.bib:
	wget -r -q -O ecological_models.bib http://127.0.0.1:23119/better-bibtex/export/collection?/2/F6ZJH5MK.biblatex

sync: sync.remember

sync.remember:
	cd ..; ./sync_remember.sh; cd overleaf;

%.tex: %.Rtex ecological_models.bib resources/* figures/*
	Rscript -e "library(knitr); knit('$<')"

%.pdf: %.tex
	latexmk -f -pdf $<

clean: 
	latexmk -C *.tex
	rm -f *.tex
	rm -f *.bbl
	rm -f *.run.xml

clean_cache:
	rm -f figures/knitr-*

viewpdf: all
	evince *.pdf

spell:
	aspell -c -t --tex-check-comments -b text.tex

pdf: all

.PHONY: clean clean_cache all ecological_models.bib sync.remember
.PRECIOUS: %.tex
