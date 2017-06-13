# this will load qtl2geno, qtl2scan, qtl2plot
library(qtl)

# convert data to R/qtl2 format
sug2 <- convert2cross2(sug)

# summary
summary(sug2)

# other summary functions
n_ind(sug2)
n_pheno(sug2)
n_mar(sug2)
tot_mar(sug2)
