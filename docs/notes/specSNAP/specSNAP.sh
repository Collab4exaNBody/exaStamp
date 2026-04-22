#!/bin/bash
rm -f *.ilg *.lof *.lot *.bbl *.blg *.idx *.ing *.ind *.log *.aux *.nav *.out *.snm *.toc *.dvi *.ps
latex     specSNAP.tex 
makeindex specSNAP.idx # Au besoin, a ajouter pour avoir l'index
bibtex    specSNAP     # Au besoin, a ajouter pour avoir la biblio
latex     specSNAP.tex
latex     specSNAP.tex
dvips     specSNAP.dvi
ps2pdf    specSNAP.ps
rm -f *.ilg *.lof *.lot *.bbl *.blg *.idx *.ing *.ind *.log *.aux *.nav *.out *.snm *.toc *.dvi *.ps

