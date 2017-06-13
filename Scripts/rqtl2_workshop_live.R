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

# plot results with the previous ones
plot(out, map)
plot(out_lmm, map, col="orchid", add=TRUE)

# calculate kinship matrices with "loco" method
k_loco <- calc_kinship(pr, "loco")

# genome scan by LMM with "loco" method
out_loco <- scan1(pr, sug2$pheno, k_loco)

plot(out_loco, map, col="green3", add=TRUE)

# X chromosome 
Xcovar <- get_x_covar(sug2)
out <- scan1(pr, sug2$pheno, Xcovar=Xcovar)

# qtl effects
eff7 <- scan1coef(pr[,7], sug2$pheno[,1])
plot(eff7, map[7], columns=1:3)

# same using contrasts matrix
contr <- cbind(mu=c(1,1,1), a=c(-0.5, 0, 0.5), d=c(-0.5, 1, -0.5))
eff7ad <- scan1coef(pr[,7], sug2$pheno[,1], contrasts=contr)
plot(eff7ad, map[7], columns=2:3)

# same using contrasts matrix
eff15ad <- scan1coef(pr[,15], sug2$pheno[,1], contrasts=contr)
plot(eff15ad, map[15], columns=2:3)

## diversity outbreds
DOex <- read_cross2("http://rqtl.org/DOex.zip")

# calc genotype probabilities
map <- insert_pseudomarkers(DOex$gmap, step=1)
pr <- calc_genoprob(DOex, map, cores=0)

# convert genotypes probs -> allele "dosages"
apr <- genoprob_to_alleleprob(pr)

# calculate kinship matrix by "loco" method
# (but note we're just using 3 chr here, and one of them is X)
k <- calc_kinship(apr, "loco")

# genome scan by LMM with "loco" method
out <- scan1(apr, DOex$pheno, k)

# plot the first one
plot(out, map)

# find peaks as before
find_peaks(out, map)

# full model with all 36 genotypes
out_full <- scan1(pr, DOex$pheno, k)
plot(out_full, map)

# qtl effects
eff2 <- scan1coef(apr[,"2"], DOex$pheno, k["2"])

plot_coefCC(eff2, map["2"])
legend("bottomleft", names(CCcolors), col=CCcolors, lwd=2, ncol=2,
    bg="gray90")

blup2 <- scan1blup(apr[,"2"], DOex$pheno, k["2"])
plot_coefCC(blup2, map["2"])
legend("bottomleft", names(CCcolors), col=CCcolors, lwd=2, ncol=2,
       bg="gray90")

# snp association tests
url <- "http://rqtl.org/c2_snpinfo.rds"
file <- basename(url)
download.file(url, file)
snpinfo <- readRDS(file)
head(snpinfo)

# calculate "strain distribution patterns" (SDP)
snpinfo$sdp <- calc_sdp(snpinfo[, 5:12])

# index snps
pmap <- interp_map(map, DOex$gmap, DOex$pmap)
snpinfo <- index_snps(pmap, snpinfo)
head(snpinfo)

# genotype probs -> snp probs
snp_pr <- genoprob_to_snpprob(pr, snpinfo)

# do the actual snp association analysis
out_snp <- scan1(snp_pr, DOex$pheno, k[["2"]])

plot(out_snp, snpinfo)

top_snps(out_snp, snpinfo, show_all_snps=FALSE)
