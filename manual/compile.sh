#!/bin/bash

pdflatex --shell-escape manual.tex
pdflatex --shell-escape manual.tex
bibtex manual.aux
pdflatex --shell-escape manual.tex
pdflatex --shell-escape manual.tex
