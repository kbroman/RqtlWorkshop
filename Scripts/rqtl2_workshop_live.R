# load R/qtl2
library(qtl2)

# load some data into R
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

# convert to qtl2 version
sug2 <- convert2cross2(sug)

# summary of object
summary(sug2)
pheno_names(sug2)
n_ind(sug2)
n_phe(sug2)
tot_mar(sug2)
n_mar(sug2)

# some things don't work:
## plot(sug2)

# calculate genotype probabilities
map <- insert_pseudomarkers(sug2$gmap, step=1)
pr <- calc_genoprob(sug2, map)

# genome scan
out <- scan1(pr, sug2$pheno)

# plot LOD curves for the first one
plot(out, map)
# add lod curves for 2nd phenotype
plot(out, map, lodcolumn=2, col="orchid", add=TRUE)

# find peaks
find_peaks(out, map)

# permutations
operm <- scan1perm(pr, sug2$pheno[,1], n_perm=1000)

# threshold
summary(operm)

# add threshold to previous plot
abline(h=summary(operm))

# calculate empirical kinship matrix
k <- calc_kinship(pr)

hist(k[lower.tri(k)], breaks=150)

# genome scan by LMM
out_lmm <- scan1(pr, sug2$pheno, k)
