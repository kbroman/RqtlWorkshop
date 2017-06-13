rqtl2_workshop.pdf: rqtl2_workshop.tex Figs/genoprobs.pdf
	pdflatex rqtl2_workshop

Figs/genoprobs.pdf: R/genoprobs.R
	cd R;R CMD BATCH $(<F)
