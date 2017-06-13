# load R/qtl
library(qtl)

# load some data into R
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

# summary of data
summary(sug)

# bits of that
nind(sug)
totmar(sug)
nphe(sug)
nmar(sug)
nchr(sug)
phenames(sug)

# summary plot
plot(sug)

# individual panels
plotMissing(sug)
plotMap(sug)
plotPheno(sug, 1)
plotPheno(sug, 2)

# calculate genotype probabilities
sug <- calc.genoprob(sug, step=1)

# genome scan for QTL
out.em <- scanone(sug)

# plot results
plot(out.em)

# biggest LOD on each chromosome
summary(out.em)

# QTL with LOD > 2
summary(out.em, threshold=2)

# interactive plots
library(qtlcharts)
iplotScanone(out.em, sug)
iplotScanone(out.em, sug, pxgtype="raw")

# permutation
operm <- scanone(sug, n.perm=200)

summary(operm)
plot(operm)

# add a title
plot(operm, main="permutations for BP")

# 5% and 20% thresholds
summary(operm, alpha=c(0.05, 0.2))

# 1.5-LOD support interval
lodint(out.em, chr=7)

# 2-LOD support interval
lodint(out.em, chr=7, drop=2)

# 2-LOD support interval, forced to expand to flanking markers
lodint(out.em, chr=7, drop=2, expandtomarkers=TRUE)

## Haley-Knott regression
out.hk <- scanone(sug, method="hk")

# plot results with result from interval mapping
plot(out.em, out.hk, col=c("slateblue", "orchid"), lty=1:2)

# plot the differences
plot(out.em - out.hk, ylim=c(-0.4, 0.4), ylab="LOD(EM) - LOD(HK)")

# permutation tests are way faster
operm.hk <- scanone(sug, method="hk", n.perm=1000)

## multiple imputation
# perform the imputations
sug <- sim.geno(sug, step=1, n.draws=32)

# QTL scan
out.imp <- scanone(sug, method="imp")

plot(out.em, out.imp, col=c("slateblue", "green3"), lty=1:2)

# effect plots
iplotScanone(out.hk, sug)

# find peak marker and plot pheno vs geno
max(out.hk)
(marker <- find.marker(sug, 7, 47.7))
plot.pxg(sug, marker=marker)

# plot ave phenotype for each genotype
effectplot(sug, mname1="7@47.7")

## let's look at heart weight
plotPheno(sug, 4)

# transform to normal quantiles
y <- pull.pheno(sug, 4)
z <- nqrank(y)

# histogram of transformed phenotypes
hist(z, breaks=30)

# scanone with phenotype 4
out.hw <- scanone(sug, pheno.col=4)
out.hw <- scanone(sug, pheno.col="heart_wt")

# scanone with transformed version
out.hw_tr <- scanone(sug, pheno.col=z)

# plot them
plot(out.hw, out.hw_tr)

# plot the differences
plot(out.hw - out.hw_tr, ylim=c(-0.4, 0.4))

# non-parametric interval mapping
out.np <- scanone(sug, pheno.col="heart_wt", 
                  model='np')

plot(out.np - out.hw_tr, ylim=c(-0.4,0.4))

## 2d, two-QTL scan: scantwo
# go to a more-coarse genotype probabilities
sug <- calc.genoprob(sug, step=5)

out2 <- scantwo(sug, method="hk")

plot(out2)
plot(out2, lower="fv1")

# to make sense of this, need scantwo permutations
operm2 <- scantwo(sug, method="hk", n.perm=1000)

# load pre-calculated scantwo results
load(url("http://rqtl.org/various.RData"))

summary(out2, perms=operm2, alpha=0.05)
summary(out2, perms=operm2, alpha=0.001)

# automated multi-QTL model building
(pen <- calc.penalties(operm2, alpha=0.05))
out.sq <- stepwiseqtl(sug, max.qtl=5, penalties=pen, method="hk")

plot(out.sq)

plotLodProfile(out.sq)
par(mfrow=c(1,1)) # <- get myself out of some mess

plot(out.hk, chr=c(7, 15), add=TRUE,col="purple")

# can get lod support intervals here too
lodint(out.sq, qtl.index=1)
lodint(out.sq, qtl.index=2)
